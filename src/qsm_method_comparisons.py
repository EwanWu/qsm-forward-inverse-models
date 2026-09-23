"""Controlled 3-D QSM comparisons. Original teaching implementations; no clinical/toolbox equivalence."""
from pathlib import Path
import json, numpy as np
from scipy.ndimage import binary_erosion
from scipy.sparse.linalg import LinearOperator,cg,lsmr
from scipy.sparse import coo_matrix
from scipy.fft import fftn,ifftn
B=Path(__file__).resolve().parent;O=B/'output';O.mkdir(exist_ok=True)
N=32;shape=(N,)*3;z,y,x=np.indices(shape)-N//2;K=np.meshgrid(*([np.fft.fftfreq(N)]*3),indexing='ij');kz,ky,kx=K;k2=kx*kx+ky*ky+kz*kz;r=np.sqrt(x*x+y*y+z*z);roi=r<11.5;core=binary_erosion(roi,iterations=3);shell=roi&~binary_erosion(roi,iterations=2)
rng=np.random.default_rng(20260913);eps=rng.normal(size=shape);CHECKS={};VOLUMES={}
def dipole(b=(0,0,1)):
 bx,by,bz=b;num=bx*bx*kx*kx+by*by*ky*ky+bz*bz*kz*kz
 for a,bb,w,ww in [(kx,ky,bx,by),(kx,kz,bx,bz),(ky,kz,by,bz)]:num+=2*w*ww*a*bb*((abs(a)!=.5)&(abs(bb)!=.5))
 D=np.zeros(shape);np.divide(num,k2,out=D,where=k2>0);D=1/3-D;D[0,0,0]=0;return D
D=dipole()
def A(v,d=D):return ifftn(fftn(v)*d).real
def image(v,mask=None):
 s=v[:,N//2,:].copy();return [[None if mask is not None and not mask[j,N//2,i] else round(float(s[j,i]),7) for i in range(N)] for j in range(N)]
def metrics(rec,truth,mask,observed=None):return {'rmse':float(np.sqrt(np.mean((rec[mask]-truth[mask])**2)))}
def lbv(total):
 inside=binary_erosion(roi);boundary=roi&~inside;ids=np.full(shape,-1,int);ids[inside]=np.arange(inside.sum());coord=np.array(np.where(inside)).T;rows=[];cols=[];vals=[];rhs=np.zeros(inside.sum())
 for ax in range(3):
  for sign in [-1,1]:
   ne=coord.copy();ne[:,ax]+=sign;ix=tuple(ne.T);unknown=inside[ix];rows.extend(np.arange(len(coord))[unknown]);cols.extend(ids[ix][unknown]);vals.extend(np.repeat(-1,unknown.sum()));rhs+=np.where(boundary[ix],total[ix],0)
 rows.extend(np.arange(len(coord)));cols.extend(np.arange(len(coord)));vals.extend(np.repeat(6,len(coord)));L=coo_matrix((vals,(rows,cols)),shape=(len(coord),len(coord))).tocsr();g,info=cg(L,rhs,rtol=1e-10,maxiter=1000);assert info==0
 bg=np.zeros(shape);bg[boundary]=total[boundary];bg[inside]=g;return (total-bg)*roi

def pdf(total):
 ext=~roi;nr=int(roi.sum());no=int(ext.sum())
 def mv(v):a=np.zeros(shape);a[ext]=v;return A(a)[roi]
 def rm(v):a=np.zeros(shape);a[roi]=v;return A(a)[ext]
 op=LinearOperator((nr,no),matvec=mv,rmatvec=rm);sol=lsmr(op,total[roi],damp=.001,atol=1e-7,btol=1e-7,maxiter=2400);assert sol[1] in (1,2),sol[1:3];a=np.zeros(shape);a[ext]=sol[0];return (total-A(a))*roi,{'iterations':int(sol[2]),'stop_code':int(sol[1]),'damp':.001,'normal_residual':float(sol[4])}

def resharp(total):
 kk=(x*x+y*y+z*z)<=9;kernel=np.fft.ifftshift(kk/kk.sum());h=1-fftn(kernel).real;valid=binary_erosion(roi,structure=kk[N//2-3:N//2+4,N//2-3:N//2+4,N//2-3:N//2+4]);alpha=.02
 def H(v):return ifftn(fftn(v)*h).real
 target=valid*H(roi*total);rhs=(H(valid*target)*roi)[roi]
 def mv(v):a=np.zeros(shape);a[roi]=v;return (H(valid*H(a))*roi)[roi]+alpha*alpha*v
 op=LinearOperator((roi.sum(),roi.sum()),matvec=mv);sol,info=cg(op,rhs,rtol=1e-8,maxiter=1000);assert info==0
 a=np.zeros(shape);a[roi]=sol;return a,valid,{'radius_vox':3,'alpha':alpha,'normal_relative_residual':float(np.linalg.norm(mv(sol)-rhs)/np.linalg.norm(rhs))}

background={}
for case,loc in [('central',0),('boundary',8)]:
 chi=(.3*((x+3)**2+y*y+z*z<9)+.18*((x-loc)**2+y*y+(z-3)**2<4))*roi
 outside=.9*((x-14)**2+(y-1)**2+(z+3)**2<9)+(-.5)*((x+14)**2+(y+1)**2+(z-4)**2<8)
 outside*=~roi;local=A(chi);bg=A(outside);total=local+bg+.001*eps
 # Low-order solid harmonics; a deliberately restricted model, not a stand-in for PDF.
 basis=np.stack([np.ones(shape),x,y,z,x*y,x*z,y*z,x*x-y*y,2*z*z-x*x-y*y],axis=-1)
 coef=np.linalg.lstsq(basis[shell],total[shell],rcond=None)[0];harm=(total-basis@coef)*roi
 pdfrec,pdfmeta=pdf(total);shr,smask,shrmeta=resharp(total);common=smask&core
 methodraw={'harmonic':(harm,roi,{}),'lbv':(lbv(total),roi,{}),'pdf':(pdfrec,roi,pdfmeta),'resharp':(shr,smask,shrmeta)};methods={}
 for name,(rec,mask,meta) in methodraw.items():
  rec=rec-np.mean((rec-local)[common]);m=metrics(rec,local,common);m.update({'coverage':float(mask.sum()/roi.sum()),'boundary_rmse':float(np.sqrt(np.mean((rec[shell]-local[shell])**2))) if np.all(mask[shell]) else None,**meta})
  methods[name]={'map':image(rec,mask),'error':image(rec-local,mask),'metrics':m};VOLUMES['background_'+case+'_'+name]=rec
 background[case]={'truth':image(local,roi),'total':image(total,roi),'source':image(chi,roi),'methods':methods,'common_coverage':float(common.sum()/roi.sum())}
 VOLUMES['background_'+case+'_truth']=local;VOLUMES['background_'+case+'_total']=total
 print('Background '+case,flush=True)
# Periodic finite differences and their true adjoint.
def grad(v):return np.stack([np.roll(v,-1,axis=i)-v for i in range(3)])
def adj(g):return sum(np.roll(g[i],1,axis=i)-g[i] for i in range(3))
L=sum(4*np.sin(np.pi*k)**2 for k in K)
u=rng.normal(size=shape);v=rng.normal(size=(3,)+shape);CHECKS['gradient_adjoint']=float(abs(np.sum(grad(u)*v)-np.sum(u*adj(v)))/max(1,abs(np.sum(grad(u)*v))))
# Phantom contains a small susceptibility lesion missing from one magnitude prior.
large=((x+4)/4)**2+(y/6)**2+(z/7)**2<1;lesion=(x-5)**2+y*y+(z-2)**2<4;neg=(x-3)**2+(y+4)**2+(z+4)**2<6
chi=.35*large+.2*lesion-.12*neg;truth=chi-chi.mean();field=A(truth);sigma=.004;observed=field+sigma*eps;F=fftn(observed);d2=D*D
morph_good=1+.4*large+.25*lesion+.2*neg;morph_bad=1+.4*large+.2*neg+.3*((x-4)**2+y*y+(z+6)**2<6)
weights={'tv':np.ones((3,)+shape),'morphology_aligned':np.where(abs(grad(morph_good))>.01,.15,1.),'morphology_mismatch':np.where(abs(grad(morph_bad))>.01,.15,1.)}

def tvsolve(lam,w,iterations=550):
 rho=.05;den=d2+rho*L;den[0,0,0]=1;rec=np.zeros(shape);zz=np.zeros((3,)+shape);dual=zz.copy();maxchange=0
 for it in range(iterations):
  new=ifftn((D*F+rho*fftn(adj(zz-dual)))/den).real;new-=new.mean();g=grad(new);old=zz.copy();arg=g+dual;zz=np.sign(arg)*np.maximum(abs(arg)-lam*w/rho,0);dual+=g-zz;maxchange=np.linalg.norm(new-rec)/max(np.linalg.norm(new),1e-12);rec=new
 primal=np.linalg.norm(g-zz);du=np.linalg.norm(rho*adj(zz-old));station=np.linalg.norm(A(A(rec)-observed)+rho*adj(dual));return rec,{'relative_last_change':float(maxchange),'primal':float(primal),'dual':float(du),'stationarity':float(station),'iterations':iterations}
priors={'truth':image(truth),'field':image(observed),'morphology_aligned':image(morph_good),'morphology_mismatch':image(morph_bad),'methods':{},'sigma_ppm':sigma,'sweep':{}}
lamgrid=np.geomspace(.000015,.003,10)
for name,w in weights.items():
 candidates=[]
 for lam in lamgrid:
  rec,meta=tvsolve(lam,w,iterations=400);res=float(np.sqrt(np.mean((A(rec)-observed)**2)));candidates.append((lam,res,rec,meta));
 chosen=min(candidates,key=lambda v:abs(v[1]-sigma));lam=chosen[0];rec,meta=tvsolve(lam,w,iterations=1100)
 res=float(np.sqrt(np.mean((A(rec)-observed)**2)));m=metrics(rec,truth,np.ones(shape,bool));m.update({'lambda':float(lam),'residual':res,'lesion_contrast':float(rec[lesion].mean()-rec[(r>11)&(r<13)].mean()),**meta})
 priors['methods'][name]={'map':image(rec),'error':image(rec-truth),'profile':rec[N//2+2,N//2,:].tolist(),'metrics':m};priors['sweep'][name]=[{'lambda':float(ll),'residual':rr,'rmse':float(np.sqrt(np.mean((cc-truth)**2)))} for ll,rr,cc,_ in candidates];VOLUMES['prior_'+name]=rec
 print('Prior '+name+' '+json.dumps(m),flush=True)
# Quadratic-gradient control uses its own noise-matched lambda, not the TV lambda.
grid=np.geomspace(.001,.3,80);candidates=[]
for lam in grid:
 den=d2+lam*lam*L;den[0,0,0]=1;rec=ifftn(D*F/den).real;candidates.append((lam,float(np.sqrt(np.mean((A(rec)-observed)**2))),rec))
lam,res,rec=min(candidates,key=lambda v:abs(v[1]-sigma));priors['methods']['gradient']={'map':image(rec),'error':image(rec-truth),'profile':rec[N//2+2,N//2,:].tolist(),'metrics':{'lambda':float(lam),'residual':res,'rmse':float(np.sqrt(np.mean((rec-truth)**2))),'lesion_contrast':float(rec[lesion].mean()-rec[(r>11)&(r<13)].mean())}}
priors['truth_profile']=truth[N//2+2,N//2,:].tolist();VOLUMES['prior_truth']=truth;VOLUMES['prior_observed']=observed
# Multiple orientations: equal total time, independent noise at sigma*sqrt(O).
cosmos={}
for name,angles in [('single',[(0,0)]),('limited',[(0,0),(12,0),(0,12)]),('wide',[(0,0),(60,0),(0,60)])]:
 dirs=[]
 for ax,ay in angles:
  ax=np.deg2rad(ax);ay=np.deg2rad(ay);dirs.append((np.sin(ax),np.sin(ay)*np.cos(ax),np.cos(ay)*np.cos(ax)))
 ds=[dipole(d) for d in dirs];den=sum(d*d for d in ds);num=np.zeros(shape,dtype=complex);so=sigma*np.sqrt(len(ds));rrng=np.random.default_rng(415)
 for d in ds:num+=d*fftn(A(truth,d)+so*rrng.normal(size=shape))
 obs=den>1e-10;rec=np.zeros(shape,dtype=complex);np.divide(num,den,out=rec,where=obs);rec=ifftn(rec).real
 # Report raw conditioning, no added spatial regularization; exact zeros are set to zero.
 info=den/len(ds);m={'rmse':float(np.sqrt(np.mean((rec-truth)**2))),'weak_fraction':float(np.mean(info.ravel()[1:]<.0025)),'noise_per_orientation':float(so),'orientations':angles}
 cosmos[name]={'map':image(rec),'error':image(rec-truth),'information':image(np.fft.fftshift(info)),'metrics':m};VOLUMES['cosmos_'+name]=rec
 print('COSMOS '+name+' '+json.dumps(m),flush=True)
# Algebraic checks independent of optimizers and parameter selection.
CHECKS['dipole_adjoint']=float(abs(np.sum(A(u)*chi)-np.sum(u*A(chi)))/max(1,abs(np.sum(u*A(chi)))))
CHECKS['constant_null']=float(abs(A(np.ones(shape))).max())
CHECKS['fft_roundtrip']=float(abs(ifftn(fftn(u)).real-u).max())
CHECKS['reshape_finite']=bool(all(np.all(np.isfinite(v)) for v in VOLUMES.values()))
assert max(v for v in CHECKS.values() if type(v)!=bool)<1e-10
out={'N':N,'seed':20260913,'background':background,'priors':priors,'cosmos':cosmos,'checks':CHECKS,'scope':'3-D periodic discrete dipole; central x-z slices. Synthetic teaching comparisons, not official toolbox benchmarks.'}
(O/'qsm-comparisons.json').write_text(json.dumps(out,separators=(',',':')))
np.savez_compressed(O/'qsm-comparison-volumes.npz',**VOLUMES)
(O/'qsm-checks.json').write_text(json.dumps(CHECKS,indent=2))
print('Saved comparisons',flush=True)


