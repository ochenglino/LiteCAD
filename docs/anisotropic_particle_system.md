# Anisotropic Particle System (PS) update

## Modified files and key functions

- `src/BaseModel.h`
  - Added per-vertex principal curvature/principal direction storage and getters.
  - Added `ComputePrincipalCurvaturesAndDirections()`.
- `src/BaseModel.cpp`
  - Implemented `ComputePrincipalCurvaturesAndDirections()` (quadratic fitting in vertex tangent plane).
  - Called curvature/direction computation in `ComputeScaleAndNormals()`.
- `src/Inputed_Triangle_Mesh_Reconstruction.h`
  - Added `projectedTriangleIds` workspace buffer for per-sample projection triangle ID.
- `src/Inputed_Triangle_Mesh_Reconstruction.cpp`
  - Updated PS energy/gradient in `evaluate()` from isotropic kernel to anisotropic `M_{ij}/Q_{ij}` kernel.
  - Reused `Distance_OBB::QueryResult::triangleID` to map sample `i` to input mesh local curvature frame.

## Formula mapping

For each sample point `i` and neighbor `j`:

- `d = x_i - x_j`
- `M_ij = M_i = [v_min, v_max, n] diag(s1^2, s2^2, 0) [v_min, v_max, n]^T`
- `Q_ij = Q_i = [v_min, v_max, n] diag(|s1|, |s2|, 0) [v_min, v_max, n]^T`

Implemented PS terms:

- Energy:
  - `E_PS += lambda_PS * exp( - (d^T M_i d) / (4 sigma^2) )`
- Gradient:
  - `∇_{x_i} E_PS += lambda_PS * (Q_i d) / (2 sigma^2) * exp( - (d^T M_i d) / (4 sigma^2) )`

Code uses the equivalent expanded forms:

- `d^T M_i d = s1^2 (d·v_min)^2 + s2^2 (d·v_max)^2`
- `Q_i d = |s1|(d·v_min)v_min + |s2|(d·v_max)v_max`

## Mapping strategy: sample `i` -> mesh curvature frame

1. In projection stage (`Distance_OBB::Query`), record closest triangle ID for each sample.
2. In PS loop, use that triangle’s 3 vertices and select the nearest vertex to `x_i`.
3. Use that vertex’s principal curvature/principal directions from `BaseModel`.
4. Fallback: if triangle ID is invalid, use global nearest vertex (`GetVertexID`).

This keeps changes local and reuses existing OBB query output.

## Negative curvature handling for sqrt

Requirement uses `s1 = sqrt(Kmin)`, `s2 = sqrt(Kmax)`.

To avoid NaN and keep computation stable:

- If `K < 0`, use `s = 0`.
- Otherwise use `s = sqrt(K)`.
- `Q` uses `|s|` as required.

## Frame normalization/orthogonalization

Before building `M_i/Q_i`, the code enforces:

- `n`, `v_min`, `v_max` normalization.
- Gram-Schmidt correction:
  - remove normal component from `v_min/v_max`,
  - remove `v_min` component from `v_max`,
  - normalize again.
- Fallback tangent-frame construction from `n` when degenerate.

## Parameter meaning

- `sigma`: particle interaction scale (kernel width), computed from surface area and sample count.
- `lambda_PS`: PS energy weight.
- Neighbor set: existing `closetPoints` from KD-tree radius search (`25 * sigma^2`).
- `M_ij` choice: always use the local frame/scales of sample point `i` (`M_{ij}=M_i`).
