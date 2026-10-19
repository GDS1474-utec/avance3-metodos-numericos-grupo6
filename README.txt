AVANCE 3 - MÉTODOS NUMÉRICOS - GRUPO 6

Ejecución en MATLAB Online:
1. Abrir codigo/Avance3_Farmacocinetica.mlx en MATLAB.
2. Seleccionar la carpeta del proyecto o la carpeta codigo como carpeta actual.
3. Presionar Run y revisar los resultados en datos/ y figuras/.

Se comparan Euler explícito, Euler implícito, RK4 y NSFD para
C0=10, Vmax=2, Km=1, T=12 con pasos 1, 0.5, 0.25, 0.125 y 0.0625.
Son parámetros normalizados elegidos para comparar los métodos.

El archivo datos/tabla_errores_matlab.csv es el resultado obtenido con MATLAB Online.
Sus veinte filas se compararon con datos/tabla_errores_referencia_python.csv.
La diferencia máxima en los errores y AUC fue inferior a 1e-14, y los conteos
coincidieron exactamente. Los tiempos de ejecución solo proceden de MATLAB.

Las figuras incluidas en el PDF se generaron en la verificación de referencia;
el código de MATLAB genera las figuras correspondientes al ejecutarlo.
Las plantillas de clase se adjuntan en plantillas_originales sin cambios.
