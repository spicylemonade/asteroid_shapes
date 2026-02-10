# Mathematical Formulation of the Lightcurve Inversion Problem

## 1. Shape Representation

### 1.1 Convex Shape: Spherical Harmonics Parameterization

For convex shape inversion, the asteroid surface is described by a radial function $r(\theta, \phi)$ expanded in real spherical harmonics \cite{kaasalainen2001a}:

$$r(\theta, \phi) = \sum_{l=0}^{l_{\max}} \sum_{m=-l}^{l} c_{lm} Y_l^m(\theta, \phi)$$

where:
- $\theta \in [0, \pi]$ is the colatitude, $\phi \in [0, 2\pi)$ is the longitude
- $Y_l^m(\theta, \phi)$ are the real spherical harmonics of degree $l$ and order $m$
- $c_{lm}$ are the shape coefficients to be optimized
- $l_{\max} = 8$ provides sufficient resolution (81 shape parameters)

The real spherical harmonics are defined as:

$$Y_l^m(\theta, \phi) = \begin{cases}
\sqrt{2} N_{lm} P_l^m(\cos\theta) \cos(m\phi) & m > 0 \\
N_{l0} P_l^0(\cos\theta) & m = 0 \\
\sqrt{2} N_{l|m|} P_l^{|m|}(\cos\theta) \sin(|m|\phi) & m < 0
\end{cases}$$

where $P_l^m$ are the associated Legendre polynomials and $N_{lm}$ is the normalization factor:

$$N_{lm} = \sqrt{\frac{(2l+1)}{4\pi} \frac{(l-m)!}{(l+m)!}}$$

The surface is triangulated on an icosphere mesh with $N_{\text{facets}}$ triangular facets. For each facet $k$, the centroid direction $(\theta_k, \phi_k)$ gives the radius $r_k = r(\theta_k, \phi_k)$, and the vertex positions are computed in Cartesian coordinates.

### 1.2 Non-Convex Shape: Vertex-Based Mesh

For non-convex (genetic algorithm) inversion, the shape is represented by a triangulated mesh with $N_v$ vertices \cite{bartczak2018}:

$$\mathbf{v}_i = (x_i, y_i, z_i) \quad i = 1, \ldots, N_v$$

Vertices are initialized on an icosphere and allowed to move radially:

$$\mathbf{v}_i = r_i \hat{\mathbf{n}}_i$$

where $r_i$ is the radial distance and $\hat{\mathbf{n}}_i$ is the initial unit direction vector. This parameterization allows concavities (craters, bifurcations) that convex models cannot capture.

For each triangular facet with vertices $\mathbf{v}_a, \mathbf{v}_b, \mathbf{v}_c$:

- **Normal vector**: $\hat{\mathbf{n}}_k = \frac{(\mathbf{v}_b - \mathbf{v}_a) \times (\mathbf{v}_c - \mathbf{v}_a)}{|(\mathbf{v}_b - \mathbf{v}_a) \times (\mathbf{v}_c - \mathbf{v}_a)|}$
- **Area**: $A_k = \frac{1}{2} |(\mathbf{v}_b - \mathbf{v}_a) \times (\mathbf{v}_c - \mathbf{v}_a)|$

## 2. Scattering Law

### 2.1 Lommel-Seeliger Law with Lambert Term

The bidirectional reflectance of each surface facet follows the combined Lommel-Seeliger + Lambert scattering law \cite{lommelseeliger1887, kaasalainen2001a}:

$$S(\mu, \mu_0, \alpha) = \frac{\mu_0}{\mu + \mu_0} f(\alpha) + c_L \mu_0$$

where:
- $\mu_0 = \cos(i)$ is the cosine of the incidence angle (angle between facet normal and Sun direction)
- $\mu = \cos(e)$ is the cosine of the emission angle (angle between facet normal and observer direction)
- $\alpha$ is the solar phase angle
- $f(\alpha)$ is the single-scattering phase function
- $c_L$ is the Lambert coefficient (typically $c_L \approx 0.1$)

The incidence and emission angles for facet $k$ are:

$$\mu_{0,k} = \hat{\mathbf{n}}_k \cdot \hat{\mathbf{e}}_{\odot}$$
$$\mu_k = \hat{\mathbf{n}}_k \cdot \hat{\mathbf{e}}_{\text{obs}}$$

where $\hat{\mathbf{e}}_{\odot}$ and $\hat{\mathbf{e}}_{\text{obs}}$ are the unit vectors toward the Sun and observer in the body-fixed frame.

### 2.2 Phase Function

The single-scattering phase function is approximated by an exponential-linear model:

$$f(\alpha) = A_0 \exp(-\alpha / D) + k\alpha + 1$$

where $A_0$, $D$, and $k$ are scattering parameters. For most S-type and C-type asteroids, a simplified linear phase correction suffices:

$$f(\alpha) = 1 - \beta \alpha$$

where $\beta \approx 0.01$ deg$^{-1}$ is the phase coefficient.

## 3. Viewing/Illumination Geometry from Orbital Elements

### 3.1 Keplerian Orbit Propagation

Given the orbital elements from MPCORB $(a, e, i, \Omega, \omega, M_0, t_0)$, the heliocentric position at epoch $t$ is computed via two-body propagation \cite{murray1999}:

1. **Mean anomaly**: $M(t) = M_0 + n(t - t_0)$, where $n = \sqrt{GM_\odot / a^3}$
2. **Eccentric anomaly** $E$ from Kepler's equation: $M = E - e \sin E$ (solved iteratively)
3. **True anomaly**: $\nu = 2\arctan\left(\sqrt{\frac{1+e}{1-e}} \tan\frac{E}{2}\right)$
4. **Heliocentric distance**: $r_h = a(1 - e\cos E)$
5. **Heliocentric ecliptic coordinates**: Rotation via $(\Omega, i, \omega + \nu)$

### 3.2 Observer Geometry

For observation epoch $t_j$:

- $\mathbf{r}_{\text{ast}}(t_j)$: heliocentric position of asteroid
- $\mathbf{r}_{\oplus}(t_j)$: heliocentric position of Earth (from ephemeris or approximation)
- $\mathbf{\Delta}_j = \mathbf{r}_{\text{ast}} - \mathbf{r}_{\oplus}$: geocentric vector to asteroid

**Phase angle**:
$$\alpha_j = \arccos\left(\frac{\mathbf{r}_{\text{ast}} \cdot \mathbf{\Delta}_j}{|\mathbf{r}_{\text{ast}}||\mathbf{\Delta}_j|}\right)$$

**Sun direction in body frame** (requires spin state):
$$\hat{\mathbf{e}}_{\odot,j} = R_{\text{spin}}(t_j)^{-1} R_{\text{ecl}}^{-1} \left(-\frac{\mathbf{r}_{\text{ast}}}{|\mathbf{r}_{\text{ast}}|}\right)$$

**Observer direction in body frame**:
$$\hat{\mathbf{e}}_{\text{obs},j} = R_{\text{spin}}(t_j)^{-1} R_{\text{ecl}}^{-1} \left(\frac{\mathbf{\Delta}_j}{|\mathbf{\Delta}_j|}\right)$$

### 3.3 Spin State and Body-Frame Transformation

The spin state is parameterized by:
- $(\lambda_p, \beta_p)$: ecliptic longitude and latitude of the spin pole
- $P$: sidereal rotation period
- $\phi_0$: rotational phase at reference epoch $t_0$

The spin axis direction in ecliptic coordinates:

$$\hat{\mathbf{s}} = (\cos\beta_p \cos\lambda_p, \cos\beta_p \sin\lambda_p, \sin\beta_p)$$

The rotation matrix $R_{\text{spin}}(t)$ transforms from body-fixed to ecliptic frame:

$$R_{\text{spin}}(t) = R_z(\phi(t)) \cdot R_{\text{pole}}$$

where $\phi(t) = \phi_0 + \frac{2\pi}{P}(t - t_0)$ is the rotational phase, and $R_{\text{pole}}$ orients the body $z$-axis along the spin pole.

## 4. Synthetic Brightness Computation

### 4.1 Visible and Illuminated Facet Selection

For each observation epoch $t_j$, the total reflected brightness is the sum over all facets that are both **visible** (facing the observer) and **illuminated** (facing the Sun):

$$L_{\text{model}}(t_j) = \sum_{k=1}^{N_{\text{facets}}} V_k(t_j) \cdot A_k \cdot S(\mu_k, \mu_{0,k}, \alpha_j)$$

where the visibility function is:

$$V_k(t_j) = \begin{cases} 1 & \text{if } \mu_k > 0 \text{ and } \mu_{0,k} > 0 \\ 0 & \text{otherwise} \end{cases}$$

### 4.2 Shadowing (Non-Convex Models)

For non-convex shapes, self-shadowing and mutual occultation must be checked via ray-casting. A facet $k$ is shadowed if a ray from its centroid toward the Sun intersects another facet. Similarly, it is occluded if a ray toward the observer is blocked. This is computationally expensive but essential for accurate non-convex modeling \cite{bartczak2018}.

### 4.3 Magnitude Conversion

The model brightness in magnitudes is:

$$m_{\text{model}}(t_j) = -2.5 \log_{10}\left(\frac{L_{\text{model}}(t_j)}{L_{\text{ref}}}\right) + m_{\text{offset}}$$

where $L_{\text{ref}}$ is a reference flux and $m_{\text{offset}}$ absorbs unknown calibration constants (fitted as free parameters for each lightcurve session).

## 5. Chi-Squared Objective Function

### 5.1 Dense Lightcurve Chi-Squared

For $N_{\text{lc}}$ dense lightcurve sessions, each containing $N_i$ data points:

$$\chi^2_{\text{dense}} = \sum_{i=1}^{N_{\text{lc}}} \sum_{j=1}^{N_i} \frac{\left(m_{\text{obs},ij} - m_{\text{model},ij} - \delta_i\right)^2}{\sigma_{ij}^2}$$

where:
- $m_{\text{obs},ij}$ is the observed magnitude of data point $j$ in lightcurve $i$
- $m_{\text{model},ij}$ is the model-predicted magnitude
- $\delta_i$ is the per-session magnitude offset (nuisance parameter, analytically marginalized)
- $\sigma_{ij}$ is the observational uncertainty

The per-session offset $\delta_i$ is analytically solved:

$$\delta_i = \frac{1}{N_i} \sum_{j=1}^{N_i} (m_{\text{obs},ij} - m_{\text{model},ij})$$

### 5.2 Sparse Data Chi-Squared

For sparse photometric data (calibrated absolute magnitudes) \cite{durech2009}:

$$\chi^2_{\text{sparse}} = \sum_{j=1}^{N_{\text{sparse}}} \frac{\left(m_{\text{obs},j} - m_{\text{model},j}\right)^2}{\sigma_j^2}$$

No per-session offsets are needed since sparse data are absolutely calibrated.

### 5.3 Combined Objective

$$\chi^2_{\text{total}} = w_{\text{dense}} \chi^2_{\text{dense}} + w_{\text{sparse}} \chi^2_{\text{sparse}}$$

where $w_{\text{dense}}$ and $w_{\text{sparse}}$ are relative weights (typically $w_{\text{sparse}} / w_{\text{dense}} \approx 0.1$ to account for different data qualities) \cite{durech2009}.

## 6. Regularization Terms

### 6.1 Smoothness Regularization

Penalizes high-frequency shape variations to prevent over-fitting:

$$R_{\text{smooth}} = \lambda_s \sum_{k} \sum_{k' \in \mathcal{N}(k)} |\mathbf{v}_k - \mathbf{v}_{k'}|^2$$

where $\mathcal{N}(k)$ is the set of vertices adjacent to vertex $k$, and $\lambda_s$ is the smoothness weight.

For spherical harmonics, smoothness is naturally enforced by limiting $l_{\max}$:

$$R_{\text{smooth}}^{\text{SH}} = \lambda_s \sum_{l=0}^{l_{\max}} \sum_{m=-l}^{l} l^2(l+1)^2 c_{lm}^2$$

### 6.2 Convexity Penalty

Encourages convex solutions by penalizing inward-pointing normals relative to the centroid:

$$R_{\text{convex}} = \lambda_c \sum_{k} \max\left(0, -\hat{\mathbf{n}}_k \cdot \frac{\mathbf{c}_k - \mathbf{o}}{|\mathbf{c}_k - \mathbf{o}|}\right)^2$$

where $\mathbf{c}_k$ is the facet centroid and $\mathbf{o}$ is the body center of mass.

### 6.3 Mesh Regularity Penalty (Non-Convex)

For genetic algorithm vertex-based meshes \cite{bartczak2018}:

$$R_{\text{mesh}} = \lambda_m \left[\sum_k (A_k - \bar{A})^2 + \sum_{\text{edges}} (l_e - \bar{l})^2\right]$$

where $A_k$ is facet area, $\bar{A}$ is mean area, $l_e$ is edge length, and $\bar{l}$ is mean edge length. This prevents degenerate triangles.

### 6.4 Total Objective

$$\mathcal{L} = \chi^2_{\text{total}} + R_{\text{smooth}} + R_{\text{convex}} + R_{\text{mesh}}$$

The free parameters to optimize are:
- Shape: $\{c_{lm}\}$ (convex) or $\{r_i\}$ (non-convex vertex radii)
- Spin: $(\lambda_p, \beta_p, P, \phi_0)$
- Scattering: $(c_L, \beta)$ (optionally fitted)
- Magnitude offsets: $\{\delta_i\}$ (analytically marginalized)

## References

- \cite{kaasalainen2001a}: Convex inversion shape parameterization and optimization
- \cite{kaasalainen2001b}: Complete inverse problem formulation
- \cite{bartczak2018}: Non-convex genetic algorithm approach (SAGE)
- \cite{durech2009}: Sparse + dense photometry combination
- \cite{durech2010sparse}: DAMIT database and inversion methodology
- \cite{viikinkoski2015}: ADAM multi-data fusion framework
- \cite{lommelseeliger1887}: Lommel-Seeliger scattering law
- \cite{hapke1981}: Hapke bidirectional reflectance model
- \cite{lomb1976}: Lomb periodogram for period search
- \cite{stellingwerf1978}: Phase dispersion minimization
