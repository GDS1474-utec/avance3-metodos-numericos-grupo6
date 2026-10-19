%% Modelo farmacocinético - Grupo 6
% Comparación de Euler, RK4 y diferencias finitas no estándar.
% dC/dt = -Vmax*C/(Km+C), con C(0) = C0.
clear; clc; close all;
format long g;

% Guardamos los resultados en la carpeta de trabajo de MATLAB.
% Si estamos dentro de "codigo", usamos la carpeta del proyecto.
carpetaTrabajo = pwd;
[~, carpetaActual] = fileparts(carpetaTrabajo);
if strcmpi(carpetaActual, 'codigo')
    salida = fileparts(carpetaTrabajo);
else
    salida = carpetaTrabajo;
end
carpetaFig = fullfile(salida, 'figuras');
carpetaDatos = fullfile(salida, 'datos');
if ~exist(carpetaFig,'dir'), mkdir(carpetaFig); end
if ~exist(carpetaDatos,'dir'), mkdir(carpetaDatos); end

%% 1. Parámetros y tamaños de paso
C0 = 10; Vmax = 2; Km = 1; T = 12;
pasos = [1, 0.5, 0.25, 0.125, 0.0625];
metodos = {'Euler explicito','Euler implicito','RK4','NSFD'};
tolPos = 1e-12*C0;

% Todas las mallas terminan en T = 12.
assert(all(abs(T./pasos-round(T./pasos)) < 1e-10));

%% 2. Cálculo de errores y AUC
% Referencia obtenida a partir de la ecuación integrada.
% Despejamos numéricamente y = log(C/C0) con fzero.
res = cell(0,9);
metricas = cell(0,7);
for ip = 1:numel(pasos)
    h = pasos(ip);
    t = linspace(0,T,round(T/h)+1)';
    cref = referencia_analitica(t,C0,Vmax,Km);
    AUCref = (Km*(C0-cref(end))+(C0^2-cref(end)^2)/2)/Vmax;
    for im = 1:numel(metodos)
        tic;
        [c,niter] = simular(t,C0,Vmax,Km,metodos{im},'Newton');
        seg = toc;
        err = max(abs(c-cref));
        numNeg = sum(c < -tolPos);
        auc = trapz(t,c);
        [cmax,idxMax] = max(c);
        res(end+1,:) = {h,metodos{im},err,numNeg,min(c),auc, ...
            abs(auc-AUCref),niter,seg}; %#ok<AGROW>
        if abs(h-0.25)<1e-12 || abs(h-1)<1e-12
            metricas(end+1,:) = {h,metodos{im},cmax,t(idxMax), ...
                c(end),auc,AUCref}; %#ok<AGROW>
        end
        if abs(h-0.25)<1e-12
            switch im
                case 1, cee = c;
                case 2, cei = c;
                case 3, crk = c;
                case 4, cns = c;
            end
            tt = t; cexact = cref;
        end
        if abs(h-1)<1e-12
            if im==1, cEEstress = c; tstress = t; end
            if im==3, cRKstress = c; end
        end
    end
end
TE = cell2table(res,'VariableNames',{'h','Metodo','Error_inf','Violaciones_positividad', ...
    'Concentracion_minima','AUC_trapecios','Error_AUC','Iteraciones_Nonlineales','Tiempo_segundos'});
TM = cell2table(metricas,'VariableNames',{'h','Metodo','Cmax','tmax', ...
    'C_final','AUC_trapecios','AUC_referencia'});
writetable(TE,fullfile(carpetaDatos,'tabla_errores_matlab.csv'));
writetable(TM,fullfile(carpetaDatos,'tabla_metricas_matlab.csv'));
disp(TE); disp(TM);

%% 3. Comparación de Newton y secante
% Newton vs Secante vs raiz algebraica estable de Euler implicito.
t = (0:0.25:T)';
[cN,iN] = simular(t,C0,Vmax,Km,'Euler implicito','Newton');
[cS,iS] = simular(t,C0,Vmax,Km,'Euler implicito','Secante');
[cQ,~] = simular(t,C0,Vmax,Km,'Euler implicito','Cuadratica');
val = table(max(abs(cN-cS)),max(abs(cN-cQ)),iN,iS, ...
 'VariableNames',{'MaxDif_Newton_Secante','MaxDif_Newton_Cuadratica', ...
 'Iteraciones_Newton','Iteraciones_Secante'});
writetable(val,fullfile(carpetaDatos,'validacion_matlab.csv'));
disp(val);
assert(max(abs(cN-cQ)) < 1e-7,'Fallo la validacion Euler implicito');
assert(all(cN>=-tolPos),'Euler implicito perdio positividad');

%% 4. Gráficas
f = figure('Color','w','Position',[100 100 1000 580]);
plot(tt,cexact,'k-','LineWidth',2.6); hold on;
plot(tt,cee,'--','LineWidth',1.5);
plot(tt,cei,'-.','LineWidth',1.5);
plot(tt,crk,':','LineWidth',2.1);
plot(tt,cns,'-','LineWidth',1.5);
grid on; xlabel('Tiempo normalizado'); ylabel('Concentracion normalizada');
title('Comparacion de trayectorias, h=0.25');
legend([{'Referencia'},metodos],'Location','northeast');
exportgraphics(f,fullfile(carpetaFig,'fig1_trayectorias_matlab.png'),'Resolution',240);

f = figure('Color','w','Position',[100 100 900 560]); hold on;
for im = 1:numel(metodos)
    filas = strcmp(TE.Metodo,metodos{im});
    loglog(TE.h(filas),TE.Error_inf(filas),'o-','LineWidth',1.7,'MarkerSize',7);
end
grid on; xlabel('Tamano de paso h'); ylabel('Error maximo absoluto');
title('Convergencia en norma infinito'); legend(metodos,'Location','northwest');
exportgraphics(f,fullfile(carpetaFig,'fig2_convergencia_matlab.png'),'Resolution',240);

f = figure('Color','w','Position',[100 100 960 540]);
plot(tstress,cEEstress,'o-','LineWidth',1.7); hold on;
plot(tstress,cRKstress,'s-','LineWidth',1.7);
plot(tstress,referencia_analitica(tstress,C0,Vmax,Km),'k--','LineWidth',2);
yline(0,'k:'); grid on; xlim([5,T]);
xlabel('Tiempo normalizado'); ylabel('Concentracion normalizada');
title('Malla gruesa h=1: verificar valores negativos');
legend({'Euler explicito','RK4','Referencia','C=0'},'Location','best');
exportgraphics(f,fullfile(carpetaFig,'fig3_positividad_matlab.png'),'Resolution',240);

%% 5. Variación de parámetros
casos = [5 2 1;10 2 1;15 2 1;10 1.5 1;10 2.5 1;10 2 0.5;10 2 2];
nombres = {'C0=5';'Base';'C0=15';'Vmax=1.5';'Vmax=2.5';'Km=0.5';'Km=2'};
ts = (0:0.25:T)';
sens = cell(0,8);
f=figure('Color','w','Position',[100 100 1000 600]); hold on;
for j=1:size(casos,1)
    c0i = casos(j,1); vi = casos(j,2); ki = casos(j,3);
    ce = referencia_analitica(ts,c0i,vi,ki);
    aucE = (ki*(c0i-ce(end))+(c0i^2-ce(end)^2)/2)/vi;
    [cexp,~] = simular(ts,c0i,vi,ki,'Euler explicito','Newton');
    [cr,~] = simular(ts,c0i,vi,ki,'RK4','Newton');
    [cn,~] = simular(ts,c0i,vi,ki,'NSFD','Newton');
    sens(end+1,:) = {nombres{j},c0i,vi,ki,aucE, ...
         max(abs(cexp-ce)),max(abs(cr-ce)),sum(cexp < -1e-12*c0i)}; %#ok<AGROW>
    plot(ts,ce,'LineWidth',1.5);
end
TS=cell2table(sens,'VariableNames',{'Caso','C0','Vmax','Km','AUC_analitica', ...
 'Error_Euler','Error_RK4','Negativos_Euler'});
writetable(TS,fullfile(carpetaDatos,'tabla_sensibilidad_matlab.csv'));
disp(TS);
legend(nombres,'Location','northeast'); grid on;
xlabel('Tiempo normalizado'); ylabel('Concentracion de referencia');
title('Sensibilidad un parametro a la vez');
exportgraphics(f,fullfile(carpetaFig,'fig4_sensibilidad_matlab.png'),'Resolution',240);

fprintf('\nTablas guardadas en: %s\n', carpetaDatos);
fprintf('Graficas guardadas en: %s\n', carpetaFig);
fprintf('Los tiempos de la tabla corresponden a esta ejecucion.\n');

%% Funciones locales
function [C,iterTot] = simular(t,C0,Vmax,Km,metodo,solucionador)
    assert(C0>0 && Vmax>0 && Km>0);
    C=zeros(numel(t),1); C(1)=C0;
    h=t(2)-t(1); iterTot=0;
    f=@(x)-Vmax*x/(Km+x);
    phi=Km*expm1(h*Vmax/Km)/Vmax;
    for k=1:numel(t)-1
        x=C(k);
        switch metodo
            case 'Euler explicito'
                C(k+1)=x+h*f(x);
            case 'Euler implicito'
                [C(k+1),numIter]=paso_implicito(x,h,Vmax,Km,C0,solucionador);
                iterTot=iterTot+numIter;
            case 'RK4'
                k1=f(x); k2=f(x+h*k1/2);
                k3=f(x+h*k2/2); k4=f(x+h*k3);
                C(k+1)=x+h*(k1+2*k2+2*k3+k4)/6;
            case 'NSFD'
                C(k+1)=x/(1+phi*Vmax/(Km+x));
            otherwise
                error('Metodo desconocido');
        end
        if ~isfinite(C(k+1)), error('Solucion no finita'); end
    end
end

function [x,it]=paso_implicito(c,h,v,km,c0,tipo)
    % En Euler implícito buscamos una raíz entre 0 y la concentración anterior.
    if c==0, x=0; it=0; return; end
    if c<0, error('El metodo implicito recibio un estado negativo'); end
    if strcmp(tipo,'Cuadratica')
        b=km+h*v-c; d=sqrt(b*b+4*km*c);
        if b>=0, x=2*km*c/(b+d); else, x=(-b+d)/2; end
        it=0; return;
    end
    F=@(z) z-c+h*v*z/(km+z);
    dF=@(z) 1+h*v*km/(km+z)^2;
    tolA=1e-12*c0; tolR=1e-10; maxit=50;
    inferior=0; superior=c;
    if strcmp(tipo,'Newton')
        xAnterior=c;
        for it=1:maxit
            nuevo=xAnterior-F(xAnterior)/dF(xAnterior);
            if ~isfinite(nuevo) || nuevo<inferior || nuevo>superior
                nuevo=(inferior+superior)/2;
            end
            if F(nuevo)>0, superior=nuevo; else, inferior=nuevo; end
            if abs(nuevo-xAnterior)<=tolA+tolR*abs(nuevo) && abs(F(nuevo))<=tolA
                x=nuevo; return;
            end
            xAnterior=nuevo;
        end
    elseif strcmp(tipo,'Secante')
        anterior=0; actual=c;
        for it=1:maxit
            fa=F(anterior); fb=F(actual);
            if fb==fa, nuevo=(inferior+superior)/2;
            else, nuevo=actual-fb*(actual-anterior)/(fb-fa); end
            if ~isfinite(nuevo) || nuevo<inferior || nuevo>superior
                nuevo=(inferior+superior)/2;
            end
            if F(nuevo)>0, superior=nuevo; else, inferior=nuevo; end
            if abs(nuevo-actual)<=tolA+tolR*abs(nuevo) && abs(F(nuevo))<=tolA
                x=nuevo; return;
            end
            anterior=actual; actual=nuevo;
        end
    else
        error('Solucionador no admitido');
    end
    error('No convergio Euler implicito en 50 iteraciones');
end

function c=referencia_analitica(t,c0,v,km)
    % La solución puede expresarse con la función W de Lambert.
    % Usamos la ecuación integrada y fzero.
    c=zeros(size(t));
    for j=1:numel(t)
        if t(j)==0, c(j)=c0; continue; end
        g=@(y) km*y+c0*expm1(y)+v*t(j);
        ymin=-(c0+v*t(j))/km-1;
        y=fzero(g,[ymin,0]);
        c(j)=c0*exp(y);
    end
end
