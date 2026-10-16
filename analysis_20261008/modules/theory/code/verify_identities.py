import os
"""Independent finite-n matrix identities, entropy/KL bounds and Renyi conditions."""
from pathlib import Path
import json,math,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
rng=np.random.default_rng(20261008)

def ld(S):
    if not np.all(np.isfinite(S)) or not np.allclose(S,S.T):raise ValueError('non-finite or asymmetric covariance')
    L=np.linalg.cholesky(S)
    return 2*float(np.log(np.diag(L)).sum())

def kl(P,Q,mp,mq):
    p=len(mp);v=mp-mq
    return .5*(np.trace(np.linalg.solve(Q,P))-p+v@np.linalg.solve(Q,v)+ld(Q)-ld(P))

def renyi(P,Q,mp,mq,a):
    assert a>0 and a!=1
    C=(1-a)*P+a*Q;v=mp-mq
    return (a*ld(Q)+(1-a)*ld(P)-ld(C))/(2*(a-1))+a/2*(v@np.linalg.solve(C,v))

def renyi_precision(P,Q,mp,mq,a):
    iP=np.linalg.inv(P);iQ=np.linalg.inv(Q);M=a*iP+(1-a)*iQ
    h=a*iP@mp+(1-a)*iQ@mq
    c=a*mp@iP@mp+(1-a)*mq@iQ@mq
    log_integral=-a/2*ld(P)-(1-a)/2*ld(Q)-.5*ld(M)-.5*(c-h@np.linalg.solve(M,h))
    return log_integral/(a-1)

errors={k:[] for k in ['covarianceUpdate','logdetLemma','entropy','correlationUpdate','generalizedEigenKL','RenyiTwoDerivations']}
finite=[]
for n in [7,30,1000]:
    for p in [2,5]:
        for repetition in range(10):
            X=rng.normal(size=(n,p));X=X@np.diag(np.arange(1,p+1))
            x=rng.normal(size=p)*np.arange(1,p+1);mu=X.mean(0);v=x-mu
            Y=np.vstack([X,x]);mun=Y.mean(0)
            for ddof in [0,1]:
                S=np.cov(X,rowvar=False,ddof=ddof);Sn=np.cov(Y,rowvar=False,ddof=ddof)
                alpha=(n-1)/n if ddof==1 else n/(n+1)
                beta=1/(n+1) if ddof==1 else n/(n+1)**2
                gamma=beta/alpha;q=float(v@np.linalg.solve(S,v))
                pred=alpha*S+beta*np.outer(v,v)
                le=p*math.log(alpha)+math.log1p(gamma*q)
                errors['covarianceUpdate'].append(float(np.max(np.abs(pred-Sn))))
                errors['logdetLemma'].append(abs(ld(Sn)-ld(S)-le))
                errors['entropy'].append(abs(.5*(ld(Sn)-ld(S))-.5*le))
                R=np.corrcoef(X,rowvar=False);z=v/np.sqrt(np.diag(S))
                Rp=(R+gamma*np.outer(z,z))/np.sqrt(np.outer(1+gamma*z*z,1+gamma*z*z))
                errors['correlationUpdate'].append(float(np.max(np.abs(Rp-np.corrcoef(Y,rowvar=False)))))
                L=np.linalg.cholesky(S);B=np.linalg.solve(L,Sn);B=np.linalg.solve(L,B.T).T
                vals=np.linalg.eigvalsh(B);delta=mun-mu
                spectral=.5*(sum(vals-1-np.log(vals))+delta@np.linalg.solve(S,delta))
                errors['generalizedEigenKL'].append(abs(spectral-kl(Sn,S,mun,mu)))
                if repetition==0:finite.append({'n':n,'p':p,'covarianceDenominator':'n-1' if ddof else 'n','q':q,
                    'signedEntropy':.5*le,'absoluteEntropy':.5*abs(le),'qZeroEntropy':(alpha**(-p)-1)/gamma})

# Taylor remainder holds in a reference-whitened relative norm, with exact mean term retained.
remainder=[]
for eta in [.8,.3,.1,.03,.001]:
    A=np.diag(np.array([1,-.7,.4,-.2])*eta);P=np.eye(4)+A;Q=np.eye(4)
    exact=kl(P,Q,np.zeros(4),np.zeros(4));quad=.25*np.trace(A@A)
    bound=eta*np.linalg.norm(A,'fro')**2/(6*(1-eta))
    assert abs(exact-quad)<=bound+1e-13
    remainder.append({'relativeSpectralNorm':eta,'exactCovarianceKL':float(exact),'quadratic':float(quad),
                      'absoluteRemainder':abs(float(exact-quad)),'provedBound':float(bound)})
renyi_cases=[]
for a in [.5,2.]:
    for j in range(20):
        T=rng.normal(size=(3,3));Q=T@T.T+np.eye(3)
        L=np.linalg.cholesky(Q);A=rng.normal(size=(3,3));A=.1*(A+A.T)/2
        P=L@(np.eye(3)+A)@L.T;mp=rng.normal(size=3);mq=rng.normal(size=3)
        v1=renyi(P,Q,mp,mq,a);v2=renyi_precision(P,Q,mp,mq,a)
        errors['RenyiTwoDerivations'].append(abs(v1-v2))
        assert v1>=-1e-12
        if j==0:renyi_cases.append({'alpha':a,'covarianceFormula':float(v1),'completedSquareFormula':float(v2)})

# Restricted zero-mean identity-reference Gaussian correlation models, not the actual full augmented KL.
local_renyi=[]
E0=np.array([[0,.5,-.2],[.5,0,.3],[-.2,.3,0]])
for t in [.1,.03,.01,.003,.001]:
    E=t*E0;R=np.eye(3)+E;raw=.5*np.linalg.norm(E,'fro')**2
    d2=renyi(np.eye(3),R,np.zeros(3),np.zeros(3),2)
    local_renyi.append({'t':t,'rawSimulationQuadratic':float(raw),'Renyi2_ref_to_new':float(d2),
                        'ratioToLeadingQuadratic':float(d2/raw),'equalMeans':True,'identityReference':True})
assert abs(local_renyi[-1]['ratioToLeadingQuadratic']-1)<.002

# Independent explicit Pearson and Renyi remainder checks, including non-small bounds.
pearson_bounds=[]
for j in range(200):
    Z=rng.normal(size=(30,4));R=np.corrcoef(Z,rowvar=False);z=rng.normal(size=4)*2
    g=float(rng.uniform(.0001,.2));g0=.2;M=float(np.max(np.abs(z)))
    T=np.outer(z,z)-R*(np.outer(z*z,np.ones(4))+np.outer(np.ones(4),z*z))/2
    Rp=(R+g*np.outer(z,z))/np.sqrt(np.outer(1+g*z*z,1+g*z*z))
    C=2.5*M**4+1.5*g0*M**6
    err=float(np.max(np.abs(Rp-R-g*T)));bound=C*g*g
    assert err<=bound+1e-13
    if j<5:pearson_bounds.append({'g':g,'M':M,'maxElementRemainder':err,'provedBound':bound})
renyi_bounds=[]
for r in local_renyi:
    E=r['t']*E0;ev=np.linalg.eigvalsh(E);d=float(np.max(np.abs(ev)))
    bound=(4/(3*(1-2*d))+1/(3*(1-d)))*float(sum(abs(ev)**3))
    err=abs(r['Renyi2_ref_to_new']-r['rawSimulationQuadratic'])
    assert err<=bound+1e-13
    renyi_bounds.append({'t':r['t'],'absoluteRemainder':err,'provedBound':bound})

# Concrete counterexample: append the empirical mean; correlation perturbation is exactly zero.
n=20;p=3;alpha=(n-1)/n;beta=1/(n+1);S=np.array([[2,.2,.3],[.2,1,.1],[.3,.1,1.5]])
counter={'n':n,'p':p,'addedPoint':'reference mean','q':0.,'rawCorrelationPerturbation':0.,'logDISCO':'-Infinity',
         'signedEntropy':.5*p*math.log(alpha),'absoluteEntropy':-.5*p*math.log(alpha),
         'directedKL_new_to_ref':.5*p*(alpha-1-math.log(alpha))}
assert counter['absoluteEntropy']>0 and counter['directedKL_new_to_ref']>0

# Unbiased finite-n augmented KL has a tiny initial decrease; monotone equivalence is not universal.
qturn=2*(n+1)/(n*(n+2));g=beta/alpha
def kq(q):return .5*(p*(alpha-1-math.log(alpha))+beta*q-math.log1p(g*q)+q/(n+1)**2)
assert kq(qturn)<kq(0)
counter['unbiasedKL_q_turningPoint']=qturn;counter['KL_atTurningPoint']=kq(qturn)

# Invalid domains are reported, never silently abs/floor/ridge/log1p transformed.
domain=[]
for label,Sbad in [('singular',np.diag([1.,0.])),('indefinite',np.diag([1.,-1.]))]:
    try:ld(Sbad);raise AssertionError('invalid covariance accepted')
    except np.linalg.LinAlgError:domain.append({'case':label,'status':'undefined for nonsingular-Gaussian logdet/inverse formula'})
try:renyi(3*np.eye(2),np.eye(2),np.zeros(2),np.zeros(2),2);raise AssertionError('invalid Renyi accepted')
except np.linalg.LinAlgError:domain.append({'case':'Renyi alpha2 C=2Q-P not SPD','status':'divergent integral; no finite closed-form score'})
X=rng.normal(size=(3,5));S=np.cov(X,rowvar=False)
domain.append({'case':'p=5,n=3','sampleCovarianceRank':int(np.linalg.matrix_rank(S)),
               'correlationFinite':bool(np.all(np.isfinite(np.corrcoef(X,rowvar=False)))),
               'status':'Pearson-based DISCO can be defined despite singular covariance, if all marginal variances and weight normalizer positive'})
result={'status':'PASS','seed':20261008,'finiteNMatrixCases':120,'maxErrors':{k:max(v) for k,v in errors.items()},
        'finiteNExamples':finite,'KL_TaylorRemainderChecks':remainder,'RenyiExamples':renyi_cases,
        'restrictedRenyiLocalChecks':local_renyi,'PearsonRemainderCases':200,'PearsonRemainderExamples':pearson_bounds,
        'RenyiRemainderChecks':renyi_bounds,'nonEquivalenceCounterexample':counter,'domainChecks':domain,
        'limits':'Derivations and matrix calculations, not biological entropy validation or replacement of current simulation definitions.'}
assert max(errors['entropy'])<1e-11 and max(errors['correlationUpdate'])<1e-11
assert max(errors['generalizedEigenKL'])<1e-10 and max(errors['RenyiTwoDerivations'])<1e-10
(ROOT/'results/identity_numeric_checks.json').write_text(json.dumps(result,indent=2)+'\n')
(ROOT/'logs/numeric_python.log').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','finiteNMatrixCases','maxErrors','nonEquivalenceCounterexample']},indent=2))
