#!/usr/bin/env python3
"""Strict inductive sensitivity: original clock training kernels; frozen test transforms."""
from pathlib import Path
import os,json,ast,pickle,hashlib,time,warnings,collections,sys,difflib
from concurrent.futures import ThreadPoolExecutor
os.umask(0o077)
OUT=Path(__file__).resolve().parent;PVT=OUT/'PRIVATE';OLD=Path(os.environ['DISCO_DATA_ROOT'])/'Results_final/05_organ_network/organ_clock_comparison/20261002_local_rebuild'
os.environ.setdefault('MPLCONFIGDIR',str(PVT/'mpl_cache'))
import numpy as np,pandas as pd,sklearn,scipy,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.impute import KNNImputer
from sklearn import config_context
from sklearn.linear_model import LassoCV,Lasso,ElasticNet
from sklearn.model_selection import StratifiedKFold,GridSearchCV
from sklearn.metrics import mean_absolute_error,mean_squared_error
from joblib import parallel_backend
from threadpoolctl import threadpool_limits

def log(*x):
    line=time.strftime('%H:%M:%S')+' '+' '.join(map(str,x))
    print(line,flush=True)
    with (OUT/'logs/02_train_clocks.log').open('a') as f:f.write(line+'\n')
def sha(b):return hashlib.sha256(b).hexdigest()
def dump(v,p):p.write_text(json.dumps(v,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
spec=json.loads((PVT/'raw_locked.json').read_text());meta=pd.read_csv(PVT/'meta.csv',dtype={'eid':str})
X=np.memmap(PVT/'raw_locked.bin',dtype=spec['dtype'],shape=tuple(spec['shape']),order='F',mode='r')
cols=spec['columns'];ntrain=spec['train_N'];ntest=spec['test_N'];train=np.arange(ntrain);test=np.arange(ntrain,ntrain+ntest)
assert len(meta)==len(X) and (meta['split'].iloc[train]=='train').all() and (meta['split'].iloc[test]=='test').all()
assert not set(meta.eid.iloc[train])&set(meta.eid.iloc[test])
fs=json.loads((OUT/'numeric/feature_sets.json').read_text());qc=pd.read_csv(OUT/'numeric/training_QC.csv')
(PVT/'models').mkdir(exist_ok=True)
dump({'Python':sys.version,'packages':{'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'scikit-learn':sklearn.__version__,'matplotlib':matplotlib.__version__}},OUT/'logs/Python_environment.json')

def method_imputation(method):
    done=PVT/(method+'_imputation_done.json');zpath=PVT/(method+'_imputed_selected.npy')
    if done.exists():
        record=json.loads(done.read_text());assert record['protocol_sha256']==sha((OUT/'PROTOCOL.txt').read_bytes())
        log('Reuse THIS-DIRECTORY validated imputation',method);return np.load(zpath,mmap_mode='r')
    distance_cols=qc.loc[qc[method+'_QC'],'protein'].tolist();ix=[cols.index(c) for c in distance_cols]
    needed=sorted(set(sum(fs[method].values(),[])));sidx=[distance_cols.index(c) for c in needed]
    raw=np.array(X[:,ix],copy=True)
    # R NA_real_ uses a signalling NaN bit pattern; use quiet NaN without changing missingness.
    raw[np.isnan(raw)]=np.nan
    if method=='Oh':
        mu=np.nanmean(raw[train],axis=0);sd=np.nanstd(raw[train],axis=0,ddof=0)
        assert np.isfinite(mu).all() and np.isfinite(sd).all() and (sd>0).all()
        with np.errstate(invalid='ignore'):raw=(raw-mu)/sd
        k=int(np.sqrt(ntrain));queries=np.arange(len(X))
        np.savez(PVT/'Oh_training_scaler.npz',mean=mu,sd=sd,columns=np.array(distance_cols))
    else:k=10;queries=test
    imp=KNNImputer(n_neighbors=k);imp.fit(raw[train]);assert imp.n_features_in_==len(distance_cols)
    Z=np.lib.format.open_memmap(zpath,mode='w+',dtype=np.float64,shape=(len(X),len(needed)))
    if method=='Goeminne':
        ss=json.loads((PVT/'Goeminne_train_selected.json').read_text())
        ZZ=np.memmap(PVT/'Goeminne_train_selected.bin',dtype=ss['dtype'],shape=tuple(ss['shape']),order='F',mode='r')
        assert ss['columns']==needed and ZZ.shape==(ntrain,len(needed));Z[train]=ZZ
    log('Frozen train-donor KNN',method,'queries',len(queries),'donors',ntrain,'features',len(distance_cols),'k',k)
    blocks=[queries[i:i+256] for i in range(0,len(queries),256)]
    def chunk(ii):
        with config_context(working_memory=192):v=imp.transform(raw[ii])
        observed=np.isfinite(raw[ii]);assert np.isfinite(v).all() and np.max(np.abs(v[observed]-raw[ii][observed]))<1e-12
        return ii,v[:,sidx]
    with threadpool_limits(limits=1),ThreadPoolExecutor(max_workers=4) as pool:
        for i,(ii,z) in enumerate(pool.map(chunk,blocks),start=1):
            Z[ii]=z;Z.flush();log(method,'KNN complete',sum(map(len,blocks[:i])),'/',len(queries))
    # A participant scored in isolation must equal its value scored in a batch.
    ii=int(test[0]);single=imp.transform(raw[[ii]])[:,sidx][0];batcherr=float(np.max(np.abs(single-Z[ii])))
    assert batcherr<1e-10 and np.isfinite(Z).all()
    with (PVT/(method+'_knn_imputer.pkl')).open('wb') as f:pickle.dump(imp,f)
    dump({'method':method,'train_N':ntrain,'test_N':ntest,'k':k,'distance_P':len(distance_cols),'selected_columns':needed,'donor_source':'locked_common_train_only','test_batch_invariance_max_error':batcherr,'protocol_sha256':sha((OUT/'PROTOCOL.txt').read_bytes()),'training_imputation':'official impute.knn' if method=='Goeminne' else 'fitted sklearn KNNImputer','test_imputation':'frozen common-training donors'},done)
    return np.load(zpath,mmap_mode='r')

OH=method_imputation('Oh');GO=method_imputation('Goeminne')
nbpath=Path(os.environ['DISCO_DATA_ROOT'])/'Results_final/05_organ_network/organ_clock_comparison/20261002_official_reproduction/vendor/Oh/UKB_Olink_3k_aging_model_github.ipynb'
nb=json.loads(nbpath.read_text());helpers=''.join(nb['cells'][6]['source']);cell=''.join(nb['cells'][8]['source'])
exec(compile(helpers,str(nbpath)+'#cell6','exec'));plt.show=lambda:plt.close('all')
fitcell=cell[cell.index('    #TRAIN MODEL\n'):cell.index('    #TEST AND SAVE\n')]
fitcell='\n'.join(x[4:] if x.startswith('    ') else x for x in fitcell.splitlines())+'\n'
original_fit=fitcell
fitcell=fitcell.replace("pathmodel = '../data_ukb/clocks/UKB_'+tissue+'_aging_model_240829.pkl'","pathmodel = str(PVT/'models'/('Oh_'+tissue+'.pkl'))")
assert 'LassoCV(alphas=alphas, cv=cv, tol=0.01,' in fitcell
(OUT/'logs/Oh_training_adapter.diff').write_text(''.join(difflib.unified_diff(original_fit.splitlines(True),fitcell.splitlines(True),fromfile='official_training_only',tofile='strict_path_adapter')))
(OUT/'logs/Oh_executed_fit.py').write_text(fitcell)
go_source=(OLD/'vendor/clock_utils.py').read_text();node=next(n for n in ast.parse(go_source).body if isinstance(n,ast.FunctionDef) and n.name=='elastic_grid_search')
go_fit='\n'.join(go_source.splitlines()[node.lineno-1:node.end_lineno])+'\n';exec(compile(go_fit,'official clock_utils.py#elastic_grid_search','exec'));(OUT/'logs/Goeminne_executed_fit.py').write_text(go_fit)
features={};predictions=[];diagnostics=[];warnings_records=[];coefrows=[];performance_CUTOFF=.95;tocoef=[]
for method,Z in [('Oh',OH),('Goeminne',GO)]:
    needed=json.loads((PVT/(method+'_imputation_done.json')).read_text())['selected_columns']
    for organ,candidates in fs[method].items():
        candidates=list(candidates) if method=='Oh' else sorted(candidates)
        ix=[needed.index(c) for c in candidates];xx=np.asarray(Z[:,ix]);xtrain=xx[train];ytrain=meta.age.to_numpy()[train];tissue=organ;plist=candidates
        path=PVT/'models'/f'{method}_{organ}.pkl';cv_score=None
        log('TRAIN ONLY clock',method,organ,'N',len(xtrain),'P',len(candidates))
        with warnings.catch_warnings(record=True) as ws,threadpool_limits(limits=1):
            warnings.simplefilter('always')
            if method=='Oh':
                exec(compile(fitcell,'Oh official training kernel','exec'));fitted=lasso;cv_score=float(-model.mse_path_.mean(axis=1).min())
            else:
                with parallel_backend('threading',n_jobs=4):grid=elastic_grid_search(pd.DataFrame(xtrain,columns=candidates),pd.Series(ytrain),validation_folds=5,n_jobs=4,l1_ratio_search=[0.,.2,.4,.8,1.])
                fitted=grid.best_estimator_;cv_score=float(grid.best_score_)
                with path.open('wb') as f:pickle.dump(fitted,f)
        convergence=bool(fitted.dual_gap_<=fitted.tol*np.var(ytrain))
        warnings_records.append({'method':method,'organ':organ,'warning_counts':dict(collections.Counter(type(w.message).__name__+': '+str(w.message) for w in ws)),'selected_converged':convergence,'selected_n_iter':int(fitted.n_iter_),'dual_gap':float(fitted.dual_gap_),'dual_gap_tolerance':float(fitted.tol*np.var(ytrain))})
        raw=fitted.predict(pd.DataFrame(xx,columns=candidates) if method=='Goeminne' else xx);independent=xx@fitted.coef_+fitted.intercept_;err=float(np.max(np.abs(raw-independent)));assert err<1e-9
        # Training-calibrated age gap; never fit calibration/SD on test.
        A=np.column_stack([np.ones(ntrain),ytrain]);beta=np.linalg.lstsq(A,raw[train],rcond=None)[0];gap=raw-np.column_stack([np.ones(len(raw)),meta.age])@beta;gap_sd=float(np.std(gap[train],ddof=1));assert gap_sd>0
        predictions.append(pd.DataFrame({'eid':meta.eid,'organ':organ,'method':method,'split':meta['split'],'score':raw,'gap_train_calibrated':gap,'gap_z_trainingSD':gap/gap_sd}))
        for protein,coefficient in zip(candidates,fitted.coef_):coefrows.append({'method':method,'organ':organ,'protein':protein,'coefficient':float(coefficient),'nonzero':bool(coefficient!=0)})
        features[f'{method}_{organ}']=[c for c,b in zip(candidates,fitted.coef_) if b!=0]
        for split,ii in [('train',train),('test',test)]:diagnostics.append({'method':method,'organ':organ,'split':split,'N':len(ii),'r':float(np.corrcoef(raw[ii],meta.age.to_numpy()[ii])[0,1]),'MAE':float(mean_absolute_error(meta.age.to_numpy()[ii],raw[ii])),'RMSE':float(np.sqrt(mean_squared_error(meta.age.to_numpy()[ii],raw[ii]))),'candidate_P':len(candidates),'selected_P':int(np.count_nonzero(fitted.coef_)),'alpha':float(fitted.alpha),'l1_ratio':1. if method=='Oh' else float(fitted.l1_ratio),'independent_score_max_error':err,'calibration_intercept':float(beta[0]),'calibration_age_slope':float(beta[1]),'training_gap_sd':gap_sd,'CV_negative_MAE_or_negative_min_MSE':cv_score})
        pd.DataFrame(diagnostics).to_csv(OUT/'numeric/age_prediction_diagnostics.csv',index=False);dump(warnings_records,OUT/'qa/clock_training_warnings.json')
        log('COMPLETE clock',method,organ,'selected P',len(features[f'{method}_{organ}']),'converged',convergence)
pd.concat(predictions,ignore_index=True).to_csv(PVT/'clock_predictions.csv',index=False)
pd.DataFrame(coefrows).to_csv(OUT/'numeric/clock_coefficients.csv',index=False);dump(features,OUT/'numeric/new_train_selected_proteins.json')
dump({'models':18,'all_models_new_common_train_only':True,'train_N':ntrain,'test_N':ntest,'test_participants_used_in_fitting':0,'all_selected_converged':all(z['selected_converged'] for z in warnings_records),'Oh_helper_sha256':sha(helpers.encode()),'Oh_original_training_cell_sha256':sha(cell.encode()),'Goeminne_elastic_grid_search_sha256':sha(go_fit.encode()),'original_paper_numeric_reproduction':False},OUT/'qa/clocks_complete.json')
assert all(z['selected_converged'] for z in warnings_records),'Selected model not converged; hold, do not export final panel.'
log('ALL 18 TRAIN-ONLY CLOCK MODELS COMPLETE')
