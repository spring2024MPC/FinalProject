# Spring2024_ME555-05_ModelPredictiveControl_FinalProject  
**Athre Hollakal, Ziyao Yin, Nicholas Morrison**


**Discritize using Euler's Methods:**

$$x = \hat{x} + T_S f(\hat{x}, u) $$
$$x = I \hat{x} + T_S (A_c \hat{x} + B_c u + D)$$ 
$$I \hat{x} + T_S A_c \hat{x} + T_S B_c u + T_S D$$
$$x = (I. + T_S A_c) \hat{x} + (T B_c) u + T_S D$$
$$T_S = 0.2$$
where

$$A_c = \begin{bmatrix} \frac{1000}{100*V_e} & 0 & \frac{-1000}{100} & 0 & 0 \\\ 0 & \frac{1000}{100*V_e} & 0 & 0 & \frac{-1000}{100^2} \\\ 0&0&0&1&0 \\\ 0&0&0&0&0 \\\ 0&0&0&0&0 \end{bmatrix}$$

$$B_c = \begin{bmatrix} 0  &  \frac{1000}{100} \\\ \frac{1}{100}   & 0 \\\ 0  & 0 \\\ 0   & \frac{1000*L}{I_z} \\\ -\frac{1}{V_e}  & 0 \end{bmatrix}$$


$$x = (I + T_S \begin{bmatrix} \frac{1000}{100*V_e} & 0 & \frac{-1000}{100} & 0 & 0 \\ 0 & \frac{1000}{100*V_e} & 0 & 0 & \frac{-1000}{100^2} \\ 0&0&0&1&0 \\ 0&0&0&0&0 \\ 0&0&0&0&0 \end{bmatrix} ) \hat{x} + (T_S \begin{bmatrix} 0 && \frac{1000}{100} \\ \frac{1}{100} && 0 \\ 0 && 0 \\ 0 && \frac{1000*L}{I_z} \\ -\frac{1}{V_e} && 0 \end{bmatrix}) \begin{bmatrix} U_1 \\ U_2 \end{bmatrix} + T_S D$$


$$x = \begin{bmatrix} \frac{1000*0.2}{100*V_e} + 1 & 0 & \frac{-1000 * 0.2}{100} & 0 & 0 \\ 0 & \frac{1000 * 0.2}{100 * V_e} + 1 & 0 & 0 & \frac{-1000 * 0.2}{100^2} \\ 0&0&1&1&0 \\ 0&0&0&1&0 \\ 0&0&0&0&1 \end{bmatrix}  \hat{x} + \begin{bmatrix} 0 && \frac{1000* 0.2}{100} \\ \frac{0.2}{100} && 0 \\ 0 && 0 \\ 0 && \frac{1000*L*0.2}{I_z} \\ -\frac{0.2}{V_e} && 0 \end{bmatrix} \begin{bmatrix} U_1 \\ U_2 \end{bmatrix} + 0.2 D$$