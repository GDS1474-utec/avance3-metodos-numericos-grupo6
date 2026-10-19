"""Replica independiente del caso ilustrativo; para comprobar las tablas sin MATLAB.
Requiere numpy, scipy, pandas, matplotlib. No sustituye ejecutar el .mlx en MATLAB.
La referencia analítica se calcula con scipy.special.lambertw.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import lambertw
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
FIG=ROOT/'figuras'; DAT=ROOT/'datos'
FIG.mkdir(exist_ok=True);DAT.mkdir(exist_ok=True)
C0=10.; V=2.; KM=1.; T=12.
H=[1.,.5,.25,.125,.0625]
METHODS=['Euler explícito','Euler implícito','RK4','NSFD']

def exacta(t,c0=C0,v=V,km=KM):
    t=np.asarray(t,dtype=float)
    z=(c0/km)*np.exp((c0-v*t)/km)
    return (km*lambertw(z).real).astype(float)

def auc_exacta(t,c0=C0,v=V,km=KM):
    cf=float(exacta([t],c0,v,km)[0]);return (km*(c0-cf)+.5*(c0*c0-cf*cf))/v

def paso_imp(c,h,v,km,c0,tipo='Newton'):
    if c == 0:return 0.,0
    assert c>0
    if tipo=='Cuadratica':
        b=km+h*v-c; d=np.sqrt(b*b+4*km*c)
        return float(2*km*c/(b+d) if b>=0 else (-b+d)/2),0
    f=lambda z: z-c+h*v*z/(km+z)
    df=lambda z: 1+h*v*km/(km+z)**2
    low=0.;high=c;ta=1e-12*c0;tr=1e-10
    if tipo=='Newton':
        old=c
        for k in range(1,51):
            nxt=old-f(old)/df(old)
            if not np.isfinite(nxt) or nxt<low or nxt>high:nxt=(low+high)/2
            if f(nxt)>0:high=nxt
            else:low=nxt
            if abs(nxt-old)<=ta+tr*abs(nxt) and abs(f(nxt))<=ta:return nxt,k
            old=nxt
    if tipo=='Secante':
        old=0.;now=c
        for k in range(1,51):
            fa,fb=f(old),f(now)
            nxt=now-fb*(now-old)/(fb-fa) if fa!=fb else (low+high)/2
            if not np.isfinite(nxt) or nxt<low or nxt>high:nxt=(low+high)/2
            if f(nxt)>0:high=nxt
            else:low=nxt
            if abs(nxt-now)<=ta+tr*abs(nxt) and abs(f(nxt))<=ta:return nxt,k
            old,now=now,nxt
    raise ValueError('No converge '+str((c,h,v,km,tipo)))

def simular(h,metodo,c0=C0,v=V,km=KM,tend=T,tipo='Newton'):
    n=round(tend/h)
    assert abs(n*h-tend)<1e-10
    t=np.linspace(0,tend,n+1)
    y=np.zeros(n+1);y[0]=c0;it=0
    phi=km*np.expm1(h*v/km)/v
    f=lambda z:-v*z/(km+z)
    for i in range(n):
        c=y[i]
        if metodo=='Euler explícito': y[i+1]=c+h*f(c)
        elif metodo=='Euler implícito':
            y[i+1],itn=paso_imp(c,h,v,km,c0,tipo);it+=itn
        elif metodo=='RK4':
            k1=f(c);k2=f(c+h*k1/2);k3=f(c+h*k2/2);k4=f(c+h*k3)
            y[i+1]=c+h*(k1+2*k2+2*k3+k4)/6
        elif metodo=='NSFD':y[i+1]=c/(1+phi*v/(km+c))
        else:raise ValueError(metodo)
        assert np.isfinite(y[i+1])
    return t,y,it


def main():
    data=[];met=[];by={}
    for h in H:
        t=np.linspace(0,T,round(T/h)+1)
        ex=exacta(t);auc=auc_exacta(T)
        for name in METHODS:
            tt,vals,niter=simular(h,name)
            err=float(max(abs(vals-ex)))
            aucn=float(np.trapezoid(vals,tt))
            neg=int(sum(vals < -1e-12*C0))
            data.append(dict(h=h,Metodo=name,Error_inf=err,Violaciones_positividad=neg,
                 Concentracion_minima=float(min(vals)),AUC_trapecios=aucn,
                 Error_AUC=abs(aucn-auc),Iteraciones_Nonlineales=niter))
            if h in (1.,.25):
                im=int(np.argmax(vals))
                met.append(dict(h=h,Metodo=name,Cmax=float(vals[im]),tmax=float(tt[im]),
                                C_final=float(vals[-1]),AUC_trapecios=aucn,AUC_referencia=auc))
            by[(h,name)]=vals
    df=pd.DataFrame(data);dm=pd.DataFrame(met)
    df.to_csv(DAT/'tabla_errores_referencia_python.csv',index=False,float_format='%.15g')
    dm.to_csv(DAT/'tabla_metricas_referencia_python.csv',index=False,float_format='%.15g')
    tN,cN,iN=simular(.25,'Euler implícito',tipo='Newton')
    _,cS,iS=simular(.25,'Euler implícito',tipo='Secante')
    _,cQ,_=simular(.25,'Euler implícito',tipo='Cuadratica')
    val=pd.DataFrame([{'MaxDif_Newton_Secante':float(max(abs(cN-cS))),
                       'MaxDif_Newton_Cuadratica':float(max(abs(cN-cQ))),
                       'Iteraciones_Newton':iN,'Iteraciones_Secante':iS}])
    val.to_csv(DAT/'validacion_referencia_python.csv',index=False)
    assert max(abs(cN-cQ))<1e-8
    assert max(abs(cN-cS))<1e-8
    assert (by[(1.,'Euler explícito')]<-1e-12*C0).any()
    assert all((by[(h,'NSFD')]>=0).all() for h in H)
    assert all((by[(h,'Euler implícito')]>=0).all() for h in H)
    assert abs(exacta([T])[0]-8.315218e-6)<1e-9
    assert abs(auc_exacta(T)-29.99999584237)<1e-8
    sens=[]; cases=[('C0=5',5,2,1),('Base',10,2,1),('C0=15',15,2,1),
                      ('Vmax=1.5',10,1.5,1),('Vmax=2.5',10,2.5,1),
                      ('Km=0.5',10,2,.5),('Km=2',10,2,2)]
    fig,ax=plt.subplots(figsize=(8.4,4.7),layout='constrained')
    st=np.linspace(0,T,49)
    for lab,c0,v,km in cases:
        exact=exacta(st,c0,v,km)
        _,ee,_=simular(.25,'Euler explícito',c0,v,km)
        _,rk,_=simular(.25,'RK4',c0,v,km)
        sens.append(dict(Caso=lab,C0=c0,Vmax=v,Km=km,AUC_analitica=auc_exacta(T,c0,v,km),
                         Error_Euler=max(abs(ee-exact)),Error_RK4=max(abs(rk-exact)),
                         Negativos_Euler=sum(ee<-1e-12*c0)))
        ax.plot(st,exact,label=lab,linewidth=1.5)
    ax.grid(alpha=.25); ax.set(xlabel='Tiempo normalizado',ylabel='Concentración normalizada',
                               title='Sensibilidad de la solución de referencia')
    ax.legend(ncol=2,fontsize=8);fig.savefig(FIG/'fig4_sensibilidad.png',dpi=260);plt.close(fig)
    ds=pd.DataFrame(sens);ds.to_csv(DAT/'tabla_sensibilidad_referencia_python.csv',index=False,float_format='%.15g')
    for i in [1]:
        t=np.arange(0,T+1e-10,.25);ex=exacta(t)
        fig,ax=plt.subplots(figsize=(8.4,4.7),layout='constrained')
        ax.plot(t,ex,'k',lw=2.3,label='Referencia analítica')
        for meth,style in zip(METHODS,['--','-.',':','-']):
            ax.plot(t,by[(.25,meth)],style,lw=1.6,label=meth)
        ax.grid(alpha=.25);ax.set(xlabel='Tiempo normalizado',ylabel='Concentración normalizada',
                                  title='Trayectorias calculadas (h = 0,25)')
        ax.legend(fontsize=8)
        fig.savefig(FIG/'fig1_trayectorias.png',dpi=260);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.4,4.7),layout='constrained')
    for name in METHODS:
        d=df[df['Metodo']==name].sort_values('h')
        ax.loglog(d.h,d.Error_inf,'o-',ms=5,lw=1.7,label=name)
    ax.grid(which='both',alpha=.3);ax.set(xlabel='Tamaño de paso h',ylabel='Error máximo absoluto',
                                         title='Estudio de refinamiento de malla')
    ax.legend(fontsize=8);fig.savefig(FIG/'fig2_convergencia.png',dpi=260);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.4,4.5),layout='constrained')
    tt=np.arange(0,T+1.,1.)
    ax.plot(tt,by[(1.,'Euler explícito')],'o-',lw=1.6,label='Euler explícito')
    ax.plot(tt,by[(1.,'RK4')],'s-',lw=1.6,label='RK4')
    ax.plot(tt,exacta(tt),'k--',lw=2.0,label='Referencia analítica')
    ax.axhline(0,ls=':',color='0.4');ax.set_xlim(5,12)
    ax.set(xlabel='Tiempo normalizado',ylabel='Concentración normalizada',
           title='Malla gruesa h = 1: positividad')
    ax.legend(fontsize=9);ax.grid(alpha=.3)
    fig.savefig(FIG/'fig3_positividad.png',dpi=260);plt.close(fig)
    print('C(T)=',exacta([T])[0],'AUC exacta=',auc_exacta(T))
    print(df.to_string(index=False,formatters={'Error_inf':lambda x:f'{x:.6g}','AUC_trapecios':lambda x:f'{x:.6g}'}))
    print('\nVALIDACION',val.to_string(index=False))
    print('\nSENSIBILIDAD',ds.to_string(index=False))

if __name__=='__main__':main()
