AVANCE 3 - METODOS NUMERICOS - GRUPO 6

El articulo editable esta en main.tex y el PDF en main.pdf.
El anexo del articulo contiene el enlace al repositorio del grupo.

Ejecucion:
1. Abrir codigo/Avance3_Farmacocinetica.mlx en MATLAB Online.
2. Establecer la carpeta principal del proyecto (o codigo/) como carpeta de trabajo.
3. Ejecutar todas las secciones (Run).
4. Consultar las tablas de datos/ y las figuras generadas en figuras/.

Modelo y parametros normalizados: C0=10, Vmax=2, Km=1, T=12.
Mallas: h = 1, 0.5, 0.25, 0.125 y 0.0625.
Metodos: Euler explicito, Euler implicito, RK4 y NSFD.

Datos obtenidos en MATLAB Online:
- datos/tabla_errores_matlab.csv: errores y AUC para cuatro metodos y cinco mallas.
- datos/tabla_sensibilidad_metodos_matlab.csv: siete escenarios para los cuatro metodos,
  con error a h=0.25 y conteos de valores negativos a h=1.
- datos/tabla_tiempos_matlab.csv: tiempo medio por metodo en milisegundos,
  estimado con 300 ejecuciones consecutivas y h=0.25.

Comparaciones independientes:
- datos/tabla_errores_referencia_python.csv
- datos/tabla_sensibilidad_metodos_referencia.csv
Los resultados de sensibilidad de MATLAB coinciden con la referencia a una
precision absoluta mejor que 1e-14. Los tiempos proceden de MATLAB Online;
no son valores independientes ni datos de rendimiento universales.

Las plantillas de clase se conservan en plantillas_originales sin cambios.
Las unidades y parametros utilizados en el ejemplo son normalizados.
