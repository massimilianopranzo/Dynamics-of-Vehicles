# Assignment 1
## Scope
Three main task were performed:
- Fit data for pure forces and pure self-aligning moment
$$
	F_i = D_{i} \, \sin(C_{i} \, \arctan(B_{i} \, \tilde{x}-E_{i}(B_{i} \, \tilde{x} \, \arctan(B_{i} \, \tilde{x})))) + S_{V,i}  \qquad \qquad
    M_{z0}  = -t \, F_y + M_{zr} \\
	\tilde{x} \in \{\kappa,\alpha\} + S_{H,i} \quad \text{and} \quad i \in \{x , y\}
$$

- Fit data to get the weighting functions for the combined forces
$$
	F_x = G_{xa} \, F_{x0} \qquad 
	F_y = G_{yk} \, F_{y0} + S_{Vyk} 
$$

- Compare the results with a much more simple approach using the practical slips to define the weighting functions
$$
	F_x =  \frac{\lvert\sigma_x\lvert}{\sqrt{\sigma_x^2 + \sigma_y^2}} F_{x0}  \qquad
	F_y = \frac{\lvert\sigma_y\lvert}{\sqrt{\sigma_x^2 + \sigma_y^2}}F_{y0} \\ \qquad 
    \sigma_x = \frac{\kappa}{1 + \kappa} \qquad 
    \sigma_y = \frac{\tan\left(\alpha\right)}{1+\kappa} \approx \frac{\alpha}{1+\kappa}
$$

## Raw data
The dataset used refers to the tyre Hoosier 18.0 x 6.0 - 10 R25B for which the nominal vetical load corresponds to 220 N and a nominal radius of 18 in (45.72 cm).
The raw data contains measurements on longitudinal and lateral forces. For the longitudinal force, measurements also with different side slip angles were performed while for the lateral force the longitudinal slip was 0 (i.e. pure condition). Both dataset were windowed in a range where the tyre pressure is about 70 psi. <br>
The vertical load has 5 fixed values: 220, 440, 700, 900, 1120 N. <br>
The camber angle has 5 fixed values: 0 $^\circ$, 1 $^\circ$, 2 $^\circ$, 3 $^\circ$, 4 $^\circ$. <br>
The side slip angle has 3 fixed values: 0 $^\circ$, -3 $^\circ$, -6 $^\circ$. <br>
## Setup for minimization
The minimization consists in three main steps:
- Fit the data in nominal conditions (220 N, 0 $^\circ$)
- Fit the data for variable vertical load
- Fit the data for variable camber angle

The cost function to be minimized is 
$$
	res =  \frac{\sum_{i=1}^{N}(f_{i}-x_i)^{2}}{\sum_{i=1}^{N}{x_i^2}} 
$$
The quality of the fittings is evaluated by the root mean square error RMS and the coefficient of determination $R^2$.
$$
	RMS = \sqrt{\frac{\sum_{i=1}^{N}(f_{i}-x_i)^{2}}{N}} \qquad \qquad
	R^2 = 1 - \frac{\sum_{i=1}^{N}(f_{i}-x_i)^{2}}{\sum_{i=1}^{N}(x_i-\bar{x})^2}
$$
A good fitting is obtained when $RMS \approx 0$ (or much smaller than the data values) and $R^2 \approx 1$.

## Example of fitting - Long. force
The stiffness tends to increase as the vertical load increases

## Example of fitting - Combined forces
The fitting of the weighting functions were done bth with the data provided for the longitudinal force where there was both dependency on longitudinal and side slip. Both weighting functions were plotted as a function of the longitudinal slip for fixed values of side slip angles
- $G_{xa}$ is 1 if $\alpha$ = 0 for any value of $\kappa$, while it shows symmetric decreasing trend starting from 1 at the extremes of the range of $\kappa$ and decreasing as $\kappa$ tends 0 for different values of $\alpha$.
- $G_{yk}$ is 1 if $\kappa$ = 0, while it shows a bell shape as $\kappa$ tends to increase or decrease. The curves are more and more attenuated as $\alpha$ tends to 0.

The simplified combined model shows a good overlapping with the combined forces obtained with the weighting functions.

## Conclusions
Accurate fittings for pure forces
The cornering stiffness looses the proportionality with the vertical load since the data results more clustered at the extremes of the range making the force less steep around the origin. The fitting was still good but not so realistic.
The data given for the self aligning were not uniformly distributed for the different values of vertical loads: for 220 N and 440 N the data were uniformly distributed along the alpha range, while for higher values there were clusters of data only around the origin and at the extremes of the range. Furthermore the data were pretty spread making the fitting with a high deviations from the data.
Finding a good guess and proper bounds was not a easy task indeed the bounds were not always set.
<br>


# Assignment 2
## Test case
Two test cases evaluated: the first constant forward velocity and linearly increasing steering angle, the second constant steering angle and linearly increasing forward velocity. In both case we perform a positive (left turn)

## Lateral load transfer
It is half the difference between the force on the right - left tyre.  <br>
The transfer of the longitudinal force is almost 0. <br>
The transfer on the lateral direction is positive both at the front and at the rear, because there is also an increase in the vertical force.

## Normalized axle characteristics
The normalized rear axle characteristics are given by the sum of the rear forces normalized by the static vertical load on the axle. <br>
The normalized front axle characteristics are given by the projection of the longitudinal and lateral forces on the tyre frame in the longitudinal and lateral directions of the vehicle frame. Then both of them are normalize by the static vertical load on the axle.<br>
$$
F_{xf} = F_{xfl} - \sin(\delta_{fl})F_{yfl} + F_{xfr} - \sin(\delta_{fr})F_{yfr}\\
$$
On the longitudinal axle there are the rear and front axle equivalent slip angles.  <br>
$$
\alpha_f = \delta - \frac{(v + \Omega Lf)}{u}\\
\alpha_r = - \frac{v - \Omega Lr}{u}
$$

## Handling diagram and understeering gradient
The handling diagram is done by performing the difference of the axle equivalent slips for the same value of the normalized force.<br>
$$
-\Delta \alpha = \delta - \rho L = \delta_H \tau_H- \rho L  
$$

When $-\Delta \alpha < 0 $ then the vehicle has an oversteering tendency.<br>
Two different regions are present, one in which there is a linear relationship between the hand and the norm. acc and another in which this relationship is not linear anymore. <br>
The regions has been interpolated with two different polynomial functions, the first one with degree 1, while the second with degree 2. For the second interpolation, the coefficient of the term of order 1 found in the first interpolation has been used. <br>
This may be also studied by looking at the derivative, called understeering gradient. Here the gradient is normalized, since the derivation is performed w.r.t. the normalized lateral acceleration and not the real one. Here the KUS has been computed in two different methods, the first one with a differentiation of the experimental data while the second using the differentiation of the theoretical relation.
$$
\tilde{K_{US}}=-\frac{d\Delta_{\alpha}}{da_y/g}\\
$$
Here there is the great doubt about the definitions of the understeering gradient.

## Yaw rate and body slip angle
For the second test case, the yaw rate gain and the body slip angle gain have been computed. <br>
The two quantities has been computed in two different ways, the first one (blue in the plot) with the definition of (i.e. yaw rate over steering angle at the wheel, and $\beta$ over steering angle at the wheel) while the second one uses the approximation that in the linear range $\Delta \alpha = -k_{US} a_y$. In the plot it is also highlighted with the red line the acceleration velocity happening at the same time as the acceleration chosen as linearity limit for the interpolation of the handling diagram. It could be seen that the two methods are very similar for velocities below the red line, while they start to be different for higher one. <br>
Using the method, we can find the critical forward velocity for an oversteering vehicle. This velocity is given by 
$$
u = \sqrt{\frac{g}{\tilde{k_{US}} }}
$$ 
This velocity is not reached evn tough there is a sign of divergency. 

## Variation of the physical parameters
### Variable front stiffness
Suspension front stiffness increases $\Rightarrow$ lateral load transfer at the front increases $\Rightarrow$ front cornering stiffness decreases $\Rightarrow$ front slip equivalent slip increases $\Rightarrow$ understeering gradient increases

### Variable front camber angle
$\gamma>0$ means negative camber $\Rightarrow$ look at what happens at the front left wheel <br>
The camber acts as an equivalent slip  ans shifts the force to the left in the slip force plot. We can generate the same amount of force for a smaller amount of slip $\Rightarrow$ less understeering $\Rightarrow$ KUS decreases

### Front toe angle
$\delta > 0$ means positive rotation for the right wheel while negative for the left. <br>
$$
F_{yf} = \sin(\delta_{fl})F_{xfl} + F_{yfl} + \sin(\delta_{fr})F_{xfr} + F_{yfr}\\
\sin(\delta_{fl})F_{xfl} < 0 \\
\sin(\delta_{fr})F_{xfr} > 0 \\
\lvert|\sin(\delta_{fl})F_{xfl} \rvert| > \lvert| \sin(\delta_{fr})F_{xfr} \rvert|
$$
There is a net decreasing of the force on the axle $\Rightarrow$ The axle has to slide more to generate the same lateral force $\Rightarrow$ KUS increases 

## Conclusions
Two steady state test were performed <br>
The vehicle with the nominal parameters tends to be understeering <br>
The investigated parameters can modify this behavior  