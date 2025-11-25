# Derivación paso a paso del modelo RBC (Capítulo 14, Campante et al. 2021)

Este documento resume y deriva las ecuaciones clave del modelo RBC de un sector que se usa en el capítulo 14 para generar la Tabla 14.2. Se parte de preferencias separables en consumo y ocio, tecnología con rendimientos constantes a escala y un proceso AR(1) para la productividad total de los factores (PTF).

## 1. Estructura básica
- **Hogares representativos**: eligen consumo $C_t$ y trabajo $H_t$ (horas o empleo efectivo), con capital predeterminado $K_t$ y acumulación $K_{t+1} = (1-\delta)K_t + I_t$.
- **Preferencias**: utilidad esperada
  $$ E_0 \sum_{t=0}^{\infty} \beta^t \left(\frac{C_t^{1-\sigma}-1}{1-\sigma} - \theta \frac{H_t^{1+\phi}}{1+\phi}\right). $$
- **Tecnología**: producción Cobb-Douglas
  $$ Y_t = A_t K_t^{\alpha} H_t^{1-\alpha}, $$
  con PTF $A_t$ que sigue
  $$ \ln A_t = (1-\rho_A)\ln A + \rho_A \ln A_{t-1} + \varepsilon_t, \quad \varepsilon_t \sim iid\,(0,\sigma_\varepsilon^2). $$
- **Mercados competitivos**: las empresas eligen insumos tomando precios como dados; salarios $w_t$ y renta del capital $r_t$ son iguales a los productos marginales.

## 2. Problema del hogar
El hogar elige $\{C_t, H_t, K_{t+1}\}$ maximizando utilidad sujeta a:
$$ C_t + I_t = w_t H_t + r_t K_t, \quad K_{t+1} = (1-\delta)K_t + I_t. $$

### 2.1. Lagrangiano y FOCs
Lagrangiano (multiplicadores $\lambda_t$):
$$ \mathcal{L} = E_0 \sum_t \beta^t \left[ u(C_t,H_t) + \lambda_t\big(w_t H_t + r_t K_t - C_t - K_{t+1} + (1-\delta)K_t\big)\right]. $$

FOCs (derivadas parciales):
1. **Consumo**: $$ \frac{\partial \mathcal{L}}{\partial C_t}: \quad u_C(C_t,H_t) = \lambda_t. $$
2. **Trabajo**: $$ \frac{\partial \mathcal{L}}{\partial H_t}: \quad -u_H(C_t,H_t) = \lambda_t w_t. $$
3. **Capital (Euler)**: $$ \frac{\partial \mathcal{L}}{\partial K_{t+1}}: \quad \lambda_t = \beta E_t\big[ \lambda_{t+1}(r_{t+1} + 1 - \delta) \big]. $$
4. **Restricción de recursos** se satisface con igualdad.

### 2.2. Condiciones con preferencias CRRA–Frisch
- $u_C = C_t^{-\sigma}$.
- $u_H = \theta H_t^{\phi}$.

Reemplazando en las FOCs:
1. **Oferta de trabajo** (intra-temporal):
   $$ \theta H_t^{\phi} = C_t^{-\sigma} w_t \;\Rightarrow\; H_t^{\phi} = \frac{w_t}{\theta} C_t^{-\sigma}. $$
   Intuición: el costo marginal de trabajar (desutilidad) se iguala al salario real valorado en términos de utilidad marginal del consumo.
2. **Euler de consumo-capital**:
   $$ 1 = \beta E_t\left[(1+r_{t+1}-\delta)\left(\frac{C_t}{C_{t+1}}\right)^\sigma\right]. $$
   Intuición: el hogar iguala el beneficio marginal de ahorrar (rendimiento esperado descontado) con el costo marginal (perder consumo hoy).

## 3. Empresas competitivas
Maximizan beneficios $\Pi_t = A_t K_t^{\alpha}H_t^{1-\alpha} - w_t H_t - r_t K_t$.

FOCs de demanda de factores:
1. $$ w_t = (1-\alpha) A_t K_t^{\alpha} H_t^{-\alpha} = (1-\alpha) \frac{Y_t}{H_t}. $$
2. $$ r_t = \alpha A_t K_t^{\alpha-1} H_t^{1-\alpha} = \alpha \frac{Y_t}{K_t}. $$
Intuición: bajo competencia perfecta, cada factor es remunerado por su producto marginal.

## 4. Equilibrio competitivo y recursos
- **Mercados de bienes**: usando $Y_t = C_t + I_t$ (sin gobierno ni comercio exterior).
- **Acumulación de capital**: $K_{t+1} = (1-\delta)K_t + I_t$.
- **Shock tecnológico**: $\ln A_t$ como AR(1).

El equilibrio es una secuencia $(C_t,H_t,K_{t+1},Y_t,w_t,r_t,A_t)$ que satisface FOCs de hogar y firma, recursos y leyes de movimiento.

## 5. Estado estacionario
En estado estacionario $A=1$ y las variables crecen a tasa cero (sin crecimiento determinista). De la Euler:
$$ 1 = \beta (1 + r - \delta) \;\Rightarrow\; r = \frac{1}{\beta} - 1 + \delta. $$
Del producto marginal del capital: $r = \alpha \frac{Y}{K}$.
Se obtiene la relación capital-producto
$$ \frac{K}{Y} = \frac{\alpha}{r}. $$
La condición de trabajo implica
$$ H^{\phi} = \frac{1-\alpha}{\theta} \frac{Y}{H} C^{-\sigma}. $$
Usando recursos $Y = C + \delta K$, se cierran $\{C,Y,K,H\}$ en términos de parámetros $(\alpha,\beta,\delta,\sigma,\phi,\theta)$.

## 6. Log-linealización (alrededor del estado estacionario)
Para simular, el libro usa aproximación log-lineal. Denote por sombreros las desviaciones logarítmicas: $x\_t = \ln X_t - \ln X$. Las relaciones aproximadas:
- Producción: $$ y_t = a_t + \alpha k_t + (1-\alpha) h_t. $$
- Salario: $$ w_t = a_t + \alpha k_t - \alpha h_t. $$
- Renta del capital: $$ r_t = a_t + (\alpha-1) k_t + (1-\alpha) h_t. $$
- Acumulación de capital (linealizada):
  $$ k_{t+1} = (1-\delta) k_t + \delta i_t. $$
- Recursos (en shares estacionarios $c\equiv C/Y$, $i\equiv I/Y$):
  $$ y_t = c\, c_t + i\, i_t. $$
- Oferta laboral: usando la FOC intra-temporal
  $$ \phi h_t = w_t - \sigma c_t. $$
- Euler linealizada:
  $$ c_t = E_t[c_{t+1}] - \left(\frac{1-\beta(1-\delta)}{\beta\alpha}\right)(E_t[y_{t+1}] - E_t[k_{t+1}]). $$
- Shock tecnológico: $$ a_t = \rho_A a_{t-1} + \varepsilon_t. $$

Estas ecuaciones en $(k_t, c_t, h_t, y_t, i_t, a_t)$ se combinan para resolver políticas lineales $c_t = a_1 k_t + a_2 a_t$, $h_t = b_1 k_t + b_2 a_t$, etc., como se implementa en los scripts de simulación.

## 7. Intuiciones clave
- **Choque de PTF** $a_t$ aumenta el producto marginal del capital y del trabajo: suben $w_t$ y $r_t$, elevando $H_t$ y $K_{t+1}$ vía la Euler, lo que amplifica la respuesta de $Y_t$.
- **Persistencia** $\rho_A$ controla cuánto dura el choque; mayor persistencia genera respuestas más amplias de inversión y capital.
- **Elasticidad Frisch** $1/\phi$ determina cuánta oferta de trabajo se ajusta ante variaciones del salario real; menor $\phi$ produce mayor volatilidad de horas y del producto.
- **Aversion al riesgo/intertemporal** $\sigma$ afecta la suavización del consumo: mayor $\sigma$ reduce la volatilidad relativa de $C_t$ frente a $Y_t$.
- **Depreciación y paciencia** $(\delta,\beta)$ fijan la razón capital-producto y, por tanto, la sensibilidad de $Y_t$ a $K_t$ y de la inversión a choques.

## 8. Relación con Tablas 14.1 y 14.2
- Las ecuaciones anteriores permiten calibrar $\alpha,\beta,\delta,\rho_A,\sigma_\varepsilon,\sigma,\phi,\theta$ con datos de Argentina y construir el estado estacionario $(C,Y,K,H)$.
- La log-linealización se usa para simular el modelo ante choques de PTF, generar series sintéticas de $Y,C,I,K,H$ y productividad, y luego calcular desviaciones cíclicas con HP para llenar la Tabla 14.2 de momentos simulados.
- La Tabla 14.1 se obtiene aplicando HP y las fórmulas de volatilidades y correlaciones a los datos observados.

## 9. Funciones de impulso-respuesta (FIR) a un choque tecnológico positivo y transitorio
Para un choque $\varepsilon_0>0$ en $a_t$, con $a_t = \rho_A a_{t-1} + \varepsilon_t$ y $|\rho_A|<1$, las FIR log-lineales resultan de resolver el sistema de $(y_t,c_t,i_t,k_t,h_t)$ con $a_t$ exógeno y $k_{t+1} = (1-\delta)k_t + \delta i_t$:

1. **Producción**: $y_t = a_t + \alpha k_t + (1-\alpha) h_t$ responde instantáneamente por $a_t$ (impacto 1 a 1) y adicionalmente por $h_t$ y, con rezago, por $k_t$.
2. **Horas**: de $\phi h_t = w_t - \sigma c_t$ y $w_t = a_t + \alpha k_t - \alpha h_t$, se obtiene
   $$ h_t = \frac{1}{\phi + \alpha}\Big( a_t + \alpha k_t - \sigma c_t \Big). $$
   Impacto positivo: sube el salario marginal $\Rightarrow$ aumentan horas; el efecto se atenúa cuando $c_t$ reacciona al alza (si $\sigma>0$).
3. **Consumo**: la Euler
   $$ c_t = E_t[c_{t+1}] - \underbrace{\left(\tfrac{1-\beta(1-\delta)}{\beta\alpha}\right)}_{\kappa} (E_t[y_{t+1}] - E_t[k_{t+1}]) $$
   implica que un choque que eleva $E_t[y_{t+1}]$ incrementa $c_t$ ya en el impacto; la suavización hace que la respuesta sea más persistente que la de $y_t$.
4. **Inversión y capital**: usando $y_t = c\, c_t + i\, i_t$ y $k_{t+1} = (1-\delta)k_t + \delta i_t$, la inversión salta fuertemente (para aprovechar la mayor productividad futura) y el capital se ajusta gradualmente, generando forma de "hump" (respuestas crecientes en los primeros períodos).
5. **Renta del capital y tasa real**: $r_t = a_t + (\alpha-1)k_t + (1-\alpha)h_t$ sube en el impacto por $a_t$ y más si $h_t$ crece; luego cae a medida que $k_t$ sube.

### Forma típica de las FIR (intuición económica)
- **$Y_t$**: salto inmediato (por $a_t$ y $H_t$) seguido de trayectoria en campana por acumulación de $K_t$.
- **$C_t$**: incremento moderado y persistente (suavización intertemporal), con menor volatilidad relativa que $Y_t$.
- **$I_t$**: mayor elasticidad; pico inicial pronunciado y posible sobre-impulso antes de normalizarse.
- **$K_t$**: ajuste lento y monótono al nuevo nivel, retornando cuando $a_t$ se disipa.
- **$H_t$**: aumenta en el impacto; con $\phi$ alto (oferta menos elástica) la respuesta es acotada; con $\phi$ bajo es mayor.
- **$r_t$**: sube al impacto y retrocede conforme $K_t$ crece.

Si se grafican (ejes: tiempo en el eje horizontal; desviaciones en log o porcentajes en el vertical), las curvas se ven como: $Y_t$ y $H_t$ con pico contemporáneo; $I_t$ con sobre-reacción inicial; $C_t$ con curva suave; $K_t$ creciente y luego convergente.

## 10. Una extensión breve: costo de ajuste de capital
Para acercar el modelo a datos donde la inversión no es extremadamente volátil, se puede añadir un costo cuadrático de ajuste:
$$ \text{costo}_t = \frac{\psi}{2}\left(\frac{I_t}{K_t} - \delta\right)^2 K_t. $$

### Ecuaciones principales
- **Restricción de recursos**: $Y_t = C_t + I_t + \text{costo}_t$.
- **Acumulación de capital**: $K_{t+1} = (1-\delta)K_t + I_t$ (igual que antes).
- **Euler modificada** (en términos reales):
  $$ 1 = \beta E_t\Bigg[ (1+r_{t+1}-\delta)\frac{u_C(C_t,H_t)}{u_C(C_{t+1},H_{t+1})} - \psi\left(\frac{I_{t+1}}{K_{t+1}}-\delta\right) \frac{I_{t+1}}{K_{t+1}} \frac{u_C(C_t,H_t)}{u_C(C_{t+1},H_{t+1})} \Bigg]. $$
  En log-lineal, la inversión responde menos que en el RBC básico porque el término de costo penaliza saltos en $I_t/K_t$.
- **FOC de inversión** (condición de Tobin-q):
  $$ q_t = 1 + \psi\left(\frac{I_t}{K_t}-\delta\right), \quad q_t = \beta E_t\left[ q_{t+1}(1-\delta) + (1+r_{t+1}-\delta) \right]. $$
  Donde $q_t$ es el valor sombra del capital (Tobin-q). En el límite $\psi \to 0$ recupera el modelo básico.

### Hallazgos e intuición
- Un $\psi>0$ **suaviza la inversión**: las FIR de $I_t$ y $K_t$ presentan menores picos y transiciones más graduales.
- **Consumo** absorbe más del choque contemporáneo (menor presión a invertir de inmediato), elevando su correlación con $Y_t$.
- **Renta del capital** fluctúa menos bruscamente porque la acumulación de $K_t$ es más lenta.
- En datos, este ajuste ayuda a igualar la volatilidad relativa de inversión y las correlaciones a lo observado en muchos países (incluida Argentina), y puede mejorar el calce de la Tabla 14.2 frente a la 14.1.

