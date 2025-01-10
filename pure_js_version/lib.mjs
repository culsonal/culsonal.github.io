
export const dist2 = ([ x1, y1 ], [ x2, y2 ]) => Math.sqrt((x1-x2)**2 + (y1-y2)**2);
export const dist3 = ([ x1, y1, z1 ], [ x2, y2, z2 ]) => Math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2);
export const dist = (p1, p2) => p1.length == 2 ? dist2(p1, p2) : dist3(p1, p2);
export const len2 = ([ x, y ]) => Math.sqrt(x*x + y*y);
export const len2_sq = ([ x, y ]) => x*x + y*y;
export const len3 = ([ x, y, z ]) => Math.sqrt(x*x + y*y + z*z);
export const len3_sq = ([ x, y, z ]) => x*x + y*y + z*z;
export const normalize2 = ([ x, y ]) => { const l = len2([ x, y ]); return [ x/l, y/l ]; };
export const normalize3 = ([ x, y, z ]) => { const l = len3([ x, y, z ]); return [ x/l, y/l, z/l ]; };
export const normalize = pt => pt.length == 2 ? normalize2(pt) : normalize3(pt);
export const add2 = ([ x1, y1 ], [ x2, y2 ]) => [ x1+x2, y1+y2 ];
export const sub2 = ([ x1, y1 ], [ x2, y2 ]) => [ x1-x2, y1-y2 ];
export const dot2 = ([ x1, y1 ], [ x2, y2 ]) => x1*x2 + y1*y2;
export const add3 = ([ x1, y1, z1 ], [ x2, y2, z2 ]) => [ x1+x2, y1+y2, z1+z2 ];
export const sub3 = ([ x1, y1, z1 ], [ x2, y2, z2 ]) => [ x1-x2, y1-y2, z1-z2 ];
export const dot3 = ([ x1, y1, z1 ], [ x2, y2, z2 ]) => x1*x2 + y1*y2 + z1*z2;
export const mix2 = ([ x1, y1 ], [ x2, y2 ], m) => [ x1*m+(1-m)*x2, y1*m+(1-m)*y2 ];
export const mix3 = ([ x1, y1, z1 ], [ x2, y2, z2 ], m) => [ x1*m+(1-m)*x2, y1*m+(1-m)*y2, z1*m+(1-m)*z2 ];
export const mix = (p1, p2, m) => p1.length == 2 ? mix2(p1, p2, m) : mix3(p1, p2, m);
export const cross = ([ v0, v1, v2 ], [ w0, w1, w2 ]) => [ v1*w2-v2*w1, v2*w0-v0*w2, v0*w1-v1*w0 ];
export const clamp = (v, min, max) => Math.max(Math.min(v, max), min);
export const clamp2 = ([ x, y ], min, max) => [ clamp(x, min, max), clamp(y, min, max) ];
export const orth = ([ x, y ]) => [ -y, x ];
export const eye_matr3 = () => [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
export const mat_vec_mul3 = ([m1, m2, m3], v) => [dot3(m1, v), dot3(m2, v), dot3(m3, v)];
export const scale3 = (s, [x, y, z]) => [s*x, s*y, s*z];
export const scale3x3 = (s, [m1, m2, m3]) => [scale3(s, m1), scale3(s, m2), scale3(s, m3)];
export const mat_mat_sum3x3 = (m1, m2) => [ add3(m1[0], m2[0]), add3(m1[1], m2[1]), add3(m1[2], m2[2]) ];
export const outer_prod3 = ([x1, y1, z1], [x2, y2, z2]) => [
  [x1*x2, x1*y2, x1*z2], [y1*x2, y1*y2, y1*z2], [z1*x2, z1*y2, z1*z2] ];
export const det2 = ([ a, b, c, d ]) => a*d - b*c;
export const det3 = ([[a, b, c], [d, e, f], [g, h, i]]) => a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
export const trace3 = ([[a, b, c], [d, e, f], [g, h, i]]) => a+e+i;
export const mat_mul_prod2x2 = ([ a, b, c, d ], [ e, f, g, h ]) => [ a*e+b*g, a*f+b*h,  c*e+d*g, c*f+d*h ];
export const mat_mat_prod3x3 = ([[a, b, c], [d, e, f], [g, h, i]], [[j, k, l], [m, n, o], [p, q, r]]) => [
  [ a*j+b*m+c*p, a*k+b*n+c*q, a*l+b*o+c*r ],
  [ d*j+e*m+f*p, d*k+e*n+f*q, d*l+e*o+f*r ],
  [ g*j+h*m+i*p, g*k+h*n+i*q, g*l+h*o+i*r ] ];
export const transpose3x3 = ([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]) => [
  [ m11, m21, m31 ], [ m12, m22, m32 ], [ m13, m23, m33 ] ];
export const mat_inv3x3 = ([[a, b, c], [d, e, f], [g, h, i]]) =>
  scale3x3(1/det3([[a, b, c], [d, e, f], [g, h, i]]), [
    [ e*i-f*h, c*h-b*i, b*f-c*e ],
    [ f*g-d*i, a*i-c*g, c*d-a*f ],
    [ d*h-e*g, b*g-a*h, a*e-b*d ] ]);

export const tri_middle2 = ([[ x1, y1 ], [ x2, y2 ], [ x3, y3 ]]) => [ (x1+x2+x3)/3, (y1+y2+y3)/3 ];
export const tri_middle3 = ([[ x1, y1, z1 ], [ x2, y2, z2 ], [ x3, y3, z3 ]]) => [
  (x1+x2+x3)/3, (y1+y2+y3)/3, (z1+z2+z3)/3 ];
export const tri_middle = (tri, dim) => dim == 2 ? tri_middle2(tri) : tri_middle3(tri);
export const tri_normal = ([ p1, p2, p3 ]) => normalize3(cross(sub3(p2, p1), sub3(p3, p1)));
export const tri_area = ([ p1, p2, p3 ]) => 1/2*len3(cross(sub3(p3, p1), sub3(p2, p1)));
export const tet_volume = ([ p1, p2, p3, p4 ]) => 1/6*dot3(cross(sub3(p2, p1), sub3(p3, p1)), sub3(p4, p1));
export const tet_volume_grads = ([ p1, p2, p3, p4 ]) => [
  cross(sub3(p4, p2), sub3(p3, p2)), cross(sub3(p3, p1), sub3(p4, p1)),
  cross(sub3(p4, p1), sub3(p2, p1)), cross(sub3(p2, p1), sub3(p3, p1)) ];

// via sympy
const tri_area_grads = ([ p1, p2, p3 ]) => {
  const [[p1x, p1y, p1z], [p2x, p2y, p2z], [p3x, p3y, p3z]] = [p1, p2, p3];
  const [ p21, p31 ] = [ sub3(p2, p1), sub3(p3, p1) ];
  const [[ p21x, p21y, p21z ], [ p31x, p31y, p31z ]] = [ p21, p31 ];
  const [ cp21x, cp21y, cp21z ] = cross(p21, p31);
  const cp21len = len3([ cp21x, cp21y, cp21z ]);
  const d1x = (-cp21y*(p2z - p3z) + cp21z*(p2y - p3y))/cp21len;
  const d1y = ( cp21x*(p2z - p3z) - cp21z*(p2x - p3x))/cp21len;
  const d1z = (-cp21x*(p2y - p3y) + cp21y*(p2x - p3x))/cp21len;
  const d2x = (-cp21y*p31z + cp21z*p31y)/cp21len;
  const d2y = ( cp21x*p31z - cp21z*p31x)/cp21len;
  const d2z = (-cp21x*p31y + cp21y*p31x)/cp21len;
  const d3x = ( cp21y*p21z - cp21z*p21y)/cp21len;
  const d3y = (-cp21x*p21z + cp21z*p21x)/cp21len;
  const d3z = ( cp21x*p21y - cp21y*p21x)/cp21len;
  return [[d1x, d1y, d1z], [d2x, d2y, d2z], [d3x, d3y, d3z]];
};

let __seed = 1;
export const set_seed = (seed) => __seed = seed;
export const srnd = () => {
  const x = Math.sin(__seed++) * 10000;
  return x - Math.floor(x);
};

export const rnd_colors = num_colors => {
  const colors = new Float32Array(num_colors*3);
  for (let i = 0; i < colors.length; i++) colors[i] = Math.random();
  return colors;
};

export const empty_2d_array = (dim1, dim2) => {
  const arr = Array(dim1);
  for (let i = 0; i < dim1; i++) arr[i] = Array(dim2).fill(0);
  return arr;
};

export const draw_tris_2d = (ctx, { pts, tris }, lwidth=1, color='#fff') => {
  const [ w, h ] = [ ctx.canvas.width, ctx.canvas.height ];
  ctx.clearRect(0, 0, w, h);
  ctx.lineWidth = lwidth; ctx.strokeStyle = color;
  for (let i = 0; i < tris.length/3; i++) {
    const [[ x1, y1 ], [ x2, y2 ], [ x3, y3 ]] = get_tri({ pts, tris, dim: 2 }, i)[0];
    ctx.beginPath();
    ctx.moveTo(x1*w, y1*h);
    ctx.lineTo(x2*w, y2*h);
    ctx.lineTo(x3*w, y3*h);
    ctx.closePath();
    ctx.stroke();
  }
};

export const draw_facing_tris_2d = (ctx, { pts, tris, facing }, lwidth=2, color='rgba(0, 255, 0, .8)') => {
  const [ w, h ] = [ ctx.canvas.width, ctx.canvas.height ];
  ctx.lineWidth = lwidth; ctx.strokeStyle = color;
  for (let i = 0; i < facing.length/3; i++) {
    const neigh_tri_indices = facing.subarray(i*3, i*3+3);
    const [ mx, my ] = tri_middle2(get_tri({ pts, tris, dim: 2 }, i)[0]);
    for (let j = 0; j < neigh_tri_indices.length; j++) {
      const ni = neigh_tri_indices[j];
      const [ nx, ny ] = tri_middle2(get_tri({ pts, tris, dim: 2 }, ni)[0]);
      ctx.beginPath();
      ctx.moveTo(mx*w, my*h);
      ctx.lineTo(nx*w, ny*h);
      ctx.stroke();
    }
  }
};

// Doesn't actually work for grid_dim_x != grid_dim_y. Todo fix
export const make_rect_geom_2d = ([grid_dim_x, grid_dim_y]=[10, 10], [x, y]=[0, 0], [sw, sh]=[1, 1]) => {
  const [ w, h, grid_pts ] = [ 1/(grid_dim_x-1), 1/(grid_dim_y-1), empty_2d_array(grid_dim_x, grid_dim_y) ];
  for (let i = 0; i < grid_dim_x; i++) for (let j = 0; j < grid_dim_y; j++)
    grid_pts[i][j] = [ (x+i*w*sw)*2-1, (y+j*h*sh)*2-1 ];

  const tris = [];
  for (let i = 0; i < grid_dim_x; i++) for (let j = 0; j < grid_dim_y; j++) {
    if (i != 0 && j != 0)
      tris.push([ j*grid_dim_y + i, (j-1)*grid_dim_y + i, j*grid_dim_y + (i-1) ]);

    if (i != grid_dim_x-1 && j != grid_dim_y-1)
      tris.push([ j*grid_dim_y + i, (j+1)*grid_dim_y + i, j*grid_dim_y + (i+1) ]);
  }

  let idx = 0;
  const pts = new Float32Array(grid_dim_x*grid_dim_y*2);
  for (let i = 0; i < grid_dim_x; i++) for (let j = 0; j < grid_dim_y; j++) {
    pts.set(grid_pts[i][j], idx*2); idx += 1;
  }

  return [ pts, new Uint16Array(tris.flat()) ];
};

export const make_rect_geom_3d = ([grid_dim_x, grid_dim_y]=[10, 10], [x, y]=[0, 0], [sw, sh]=[1, 1]) => {
  const [ pts_2d, tris ] = make_rect_geom_2d([ grid_dim_x, grid_dim_y ], [ x, y ], [ sw, sh ]);
  const pts = new Float32Array(pts_2d.length/2 * 3);
  for (let i = 0; i < pts_2d.length/2; i++) {
    const [ x, y ] = pts_2d.subarray(i*2, i*2+2);
    pts.set([ x, y, 0 ], i*3);
  }
  
  return [ pts, tris ];
};

export const make_cylinder_geom_3d = ([gdx, gdy]=[10, 10], [x, y]=[0, 0], [sw, sh]=[.5, .5], height=.2) => {
  const [ pts, tris ] = make_rect_geom_3d([ gdx, gdy ], [ 0, 0 ], [ .5, height ]);
  for (let i = 0; i < pts.length/3; i++) {
    const ang = pts[i*3+0]*2*3.141592*2;
    pts[i*3+0] = sw*Math.sin(ang);
    pts[i*3+2] = sh*Math.cos(ang)*2;
  }

  for (let i = 0; i < pts.length/3; i++) pts[i*3+0] += x;
  for (let i = 0; i < pts.length/3; i++) pts[i*3+1] += y;

  return make_pts_tris_via_nearby_check({ pts, tris, dim: 3 });
};

export const subdivide_and_normalize = ({ pts, tris }, radius=.3) => {
  [ pts, tris ] = subdivide({ pts, tris, dim:3 });
  for (let i = 0; i < pts.length/3; i++) pts.set(scale3(radius, normalize3(pts.subarray(i*3, i*3+3))), i*3);
  return [ pts, tris ]
};

export const make_sphere_geom_3d = (radius=.3, num_subdivides=3) => {
  if (num_subdivides >= 4) { console.error('more than 3 subdivides not recommended'); return; }
  let pts  = new Float32Array([ 0, -1, 0,   1, 0, 0,   0, 1, 0,   0,  0, -1,  -1, 0, 0,   0, 0, 1 ]);
  let tris = new Uint16Array([ 5, 0, 1,  0, 1, 3,  1, 2, 5,  1, 2, 3,  4, 0, 5,  5, 2, 4,  3, 2, 4,  0, 3, 4 ]);
  for (let itr = 0; itr < num_subdivides+1; itr++) [ pts, tris ] = subdivide_and_normalize({ pts, tris }, radius);
  return [ pts, tris ];
};

export const make_half_sphere_geom_3d = (radius=.3, num_subdivides=3) => {
  let pts  = new Float32Array([ 0, -1, 0,   1, 0, 0,   0, 1, 0,   -1, 0, 0,   0, 0, -1 ]);
  let tris = new Uint16Array([  0,  1, 4,   1, 2, 4,   0, 4, 3,    2, 3, 4 ]);
  for (let itr = 0; itr < num_subdivides+1; itr++) [ pts, tris ] = subdivide_and_normalize({ pts, tris, }, radius);
  return [ pts, tris ];
};

export const update_normals = ({ pts, tris, normals }) => {
  const normals_helper = Array(normals.length/3).fill(0);
  normals.fill(0);
  for (let i = 0; i < tris.length/3; i++) {
    const [ tri, [ i1, i2, i3 ] ] = get_tri({ pts, tris, dim: 3 }, i);
    const [ nx, ny, nz ] = tri_normal(tri);
    normals[i1*3+0] += nx; normals[i1*3+1] += ny; normals[i1*3+2] += nz;
    normals[i2*3+0] += nx; normals[i2*3+1] += ny; normals[i2*3+2] += nz;
    normals[i3*3+0] += nx; normals[i3*3+1] += ny; normals[i3*3+2] += nz;
    normals_helper[i1] += 1; normals_helper[i2] += 1; normals_helper[i3] += 1;
  }

  for (let i = 0; i < normals.length/3; i++) {
    normals[i*3+0] /= normals_helper[i];
    normals[i*3+1] /= normals_helper[i];
    normals[i*3+2] /= normals_helper[i];
  }
};

// could be sped up with with a hash (for the inner loop)
export const make_pts_tris_via_nearby_check = ({ pts, tris, dim }) => {
  const new_pts = [];

  const old_indices_to_new = Array(pts.length/dim);
  const idxs_has_been_touched = Array(pts.length/dim).fill(false);

  for (let i = 0; i < pts.length/dim; i++) {
    if (idxs_has_been_touched[i]) continue;
    const p1 = pts.subarray(i*dim, i*dim+dim);
    new_pts.push(...p1);
    const new_idx = new_pts.length/dim-1;

    idxs_has_been_touched[i] = true;
    old_indices_to_new[i] = new_idx;
    for (let j = 0; j < pts.length/dim; j++) {
      const p2 = pts.subarray(j*dim, j*dim+dim);
      const eps = .0001;
      if (dist(p1, p2) < eps) {
        idxs_has_been_touched[j] = true;
        old_indices_to_new[j] = new_idx;
      }
    }
  }

  for (let i = 0; i < tris.length; i++) tris[i] = old_indices_to_new[tris[i]];
  return [ new Float32Array(new_pts), tris ];
};

export const subdivide = ({ pts, tris, dim=3 }) => {
  const [ new_pts, new_tris ] = [[], []];
  for (let i = 0; i < tris.length/3; i++) {
    const [ [ p1, p2, p3 ], [ i1, i2, i3 ] ] = get_tri({ pts, tris, dim }, i);
    const [ mid1, mid2, mid3 ] = [ mix(p1, p2, .5), mix(p2, p3, .5), mix(p3, p1, .5) ];
    new_pts.push(...p1); const p1i = new_pts.length/3-1; new_pts.push(...mid1); const m1i = new_pts.length/3-1;
    new_pts.push(...p2); const p2i = new_pts.length/3-1; new_pts.push(...mid2); const m2i = new_pts.length/3-1;
    new_pts.push(...p3); const p3i = new_pts.length/3-1; new_pts.push(...mid3); const m3i = new_pts.length/3-1;
    new_tris.push(...[ m1i, m2i, m3i,   p1i, m1i, m3i,   m1i, p2i, m2i,   m3i, m2i, p3i ]);
  }

  return make_pts_tris_via_nearby_check({ pts: new Float32Array(new_pts), tris: new Uint16Array(new_tris), dim });
};

export const make_normals = ({ pts, tris }) => {
  const normals = new Float32Array(pts.length);
  update_normals({ pts, tris, normals });
  return normals;
};

export const compute_facing_from_tris = tris => {
  const facing = new Int32Array(tris.length).fill(-1);
  const tris_to_update = Array(tris.length/3);
  for (let i = 0; i < tris_to_update.length; i++) tris_to_update[i] = i;
  update_facing_data({ tris, facing }, tris_to_update);
  return facing;
};

export const update_facing_data = ({ tris, facing }, tris_to_update) => {
  for (let i = 0; i < tris_to_update.length; i++) {
    let fidx = 0;
    const ti = tris_to_update[i];
    const tri = tris.subarray(ti*3, ti*3+3);
    for (let j = 0; j < tris_to_update.length; j++) {
      if (i == j) continue;
      const tj = tris_to_update[j];
      const tri_neigh = tris.subarray(tj*3, tj*3+3);
      const shared_face = does_tris_share_exactly_one_face(tri, tri_neigh);
      if (!shared_face) continue;
      facing[ti*3+fidx] = tj;
      fidx += 1;
    }

    // if there's any remaining facing elements that point out of our list then we want to
    // keep them; if any of them is in our list then discard since we've updated that data

    for (; fidx < 3; fidx++) for (let j = 0; j < tris_to_update.length; j++)
      if (facing[ti*3+fidx] == tris_to_update[j]) facing[ti*3+fidx] = -1;
  }
};

export const does_tris_share_exactly_one_face = (tA, tB) => {
  let sum_same_verts = 0;
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) if (tA[i] == tB[j]) sum_same_verts++;
  if (sum_same_verts != 2) return false;
  else return true;
};

export const get_tri_with_orig = ({ pts, orig, tris, dim=2 }, i) => {
  const [ i1, i2, i3 ] = tris.subarray(i*3, i*3+3);
  return [
    [  pts.subarray(i1*dim, i1*dim+dim),  pts.subarray(i2*dim, i2*dim+dim),  pts.subarray(i3*dim, i3*dim+dim) ],
    [ orig.subarray(i1*dim, i1*dim+dim), orig.subarray(i2*dim, i2*dim+dim), orig.subarray(i3*dim, i3*dim+dim) ],
    [ i1, i2, i3 ]
  ];
};

export const get_tri = ({ pts, tris, dim=2 }, i) => {
  const [ i1, i2, i3 ] = tris.subarray(i*3, i*3+3);
  return [
    [  pts.subarray(i1*dim, i1*dim+dim),  pts.subarray(i2*dim, i2*dim+dim),  pts.subarray(i3*dim, i3*dim+dim) ],
    [ i1, i2, i3 ]
  ];
};

const compute_neohookean_2d_deltas_from_tri_and_orig = (tri_curr, tri_orig, [ ch_damping, cd_damping ]) => {
  const [[ x1,  y1], [ x2,  y2], [ x3,  y3]] = tri_curr;
  const [[ox1, oy1], [ox2, oy2], [ox3, oy3]] = tri_orig;

  const X    = [   x2-x1,   x3-x1,   y2-y1,   y3-y1 ];
  const Xhat = [ ox2-ox1, ox3-ox1, oy2-oy1, oy3-oy1 ];
  const detB = 1/det2(Xhat);
  const [ b1, b2, b3, b4 ] = [ detB*Xhat[3], -detB*Xhat[1], -detB*Xhat[2], detB*Xhat[0] ];
  const [ f1, f2, f3, f4 ] = mat_mul_prod2x2(X, [ b1, b2, b3, b4 ]);
  const [ dCHdX1x, dCHdX1y ] = [ detB*(y2 - y3), detB*(x3 - x2) ];
  const [ dCHdX2x, dCHdX2y ] = [ detB*(y3 - y1), detB*(x1 - x3) ];
  const [ dCHdX3x, dCHdX3y ] = [ detB*(y1 - y2), detB*(x2 - x1) ];
  const ch_grad_sum_sq = dCHdX1x**2 + dCHdX1y**2 + dCHdX2x**2 + dCHdX2y**2 + dCHdX3x**2 + dCHdX3y**2;
  const ch_constraint = det2([ f1, f2, f3, f4 ]) - (1+cd_damping/ch_damping);
  const ch_pbd_lambda = -ch_constraint*ch_damping/ch_grad_sum_sq;
  const [ dx1h, dy1h ] = [ ch_pbd_lambda*dCHdX1x, ch_pbd_lambda*dCHdX1y ];
  const [ dx2h, dy2h ] = [ ch_pbd_lambda*dCHdX2x, ch_pbd_lambda*dCHdX2y ];
  const [ dx3h, dy3h ] = [ ch_pbd_lambda*dCHdX3x, ch_pbd_lambda*dCHdX3y ];

  const cd_constraint = Math.sqrt(f1**2 + f2**2 + f3**2 + f4**2);
  const dCDdX1x = -(f1*b1 + f1*b3 + f2*b2 + f2*b4)/cd_constraint;
  const dCDdX1y = -(f3*b1 + f3*b3 + f4*b2 + f4*b4)/cd_constraint;
  const [ dCDdX2x, dCDdX2y ] = [ (f1*b1 + f2*b2)/cd_constraint, (f3*b1 + f4*b2)/cd_constraint ];
  const [ dCDdX3x, dCDdX3y ] = [ (f1*b3 + f2*b4)/cd_constraint, (f3*b3 + f4*b4)/cd_constraint ];
  const cd_grad_sum_sq = dCDdX1x**2 + dCDdX1y**2 + dCDdX2x**2 + dCDdX2y**2 + dCDdX3x**2 + dCDdX3y**2;
  const cd_pbd_lambda = -cd_constraint*cd_damping/cd_grad_sum_sq;
  const [ dx1d, dy1d ] = [ cd_pbd_lambda*dCDdX1x, cd_pbd_lambda*dCDdX1y ];
  const [ dx2d, dy2d ] = [ cd_pbd_lambda*dCDdX2x, cd_pbd_lambda*dCDdX2y ];
  const [ dx3d, dy3d ] = [ cd_pbd_lambda*dCDdX3x, cd_pbd_lambda*dCDdX3y ];

  return [[ dx1h+dx1d, dy1h+dy1d ], [ dx2h+dx2d, dy2h+dy2d ], [ dx3h+dx3d, dy3h+dy3d ]];
};

export const solve_neohookean_constraints_tris_2d = ({ pts, tris, orig }, dampings) => {
  for (let i = 0; i < tris.length/3; i++) {
    const [ tri_curr, tri_orig, [ i1, i2, i3 ] ] = get_tri_with_orig({ pts, orig, tris, dim: 2 }, i);
    const [ [dx1, dy1], [dx2, dy2], [dx3, dy3]
    ] = compute_neohookean_2d_deltas_from_tri_and_orig(tri_curr, tri_orig, dampings);
    pts[i1*2+0] += dx1; pts[i1*2+1] += dy1;
    pts[i2*2+0] += dx2; pts[i2*2+1] += dy2;
    pts[i3*2+0] += dx3; pts[i3*2+1] += dy3;
  }
};

export const solve_neohookean_constraints_tris_3d_xy_plane = ({ pts, tris, orig }, dampings) => {
  for (let i = 0; i < tris.length/3; i++) {
    const [ tri_curr, tri_orig, [ i1, i2, i3 ] ] = get_tri_with_orig({ pts, orig, tris, dim: 3 }, i);
    const [ [dx1, dy1], [dx2, dy2], [dx3, dy3]
    ] = compute_neohookean_2d_deltas_from_tri_and_orig(tri_curr, tri_orig, dampings);
    pts[i1*3+0] += dx1; pts[i1*3+1] += dy1;
    pts[i2*3+0] += dx2; pts[i2*3+1] += dy2;
    pts[i3*3+0] += dx3; pts[i3*3+1] += dy3;
  }
};

export const solve_length_constraints_tris_2d = ({ pts, tris, orig }, damping=.1) => {
  for (let i = 0; i < tris.length/3; i++) {
    const [ curr_tri, orig_tri, tri_idxs ] = get_tri_with_orig({ pts, orig, tris, dim: 2 }, i);
    const idxs = [ [ 0, 1 ], [ 1, 2 ], [ 2, 0 ] ];
    for (let j = 0; j < idxs.length; j++) {
      const [ j1, j2 ] = idxs[j];
      const [[ x1, y1 ], [ x2, y2 ]] = [ curr_tri[j1], curr_tri[j2] ];
      const [ dist, desired_dist ] = [ dist2(curr_tri[j1], curr_tri[j2]), dist2(orig_tri[j1], orig_tri[j2]) ];
      const [ d1x, d1y ] = [ (x1-x2)/dist, (y1-y2)/dist ];
      const [ d2x, d2y ] = [ -d1x, -d1y ];
      const grad_sum_sq = d1x**2 + d1y**2 + d2x**2 + d2y**2;
      const constraint = dist - desired_dist;
      const lambda = -constraint*damping/grad_sum_sq;
      const [ i1, i2 ] = [ tri_idxs[j1], tri_idxs[j2] ];
      pts[i1*2+0] += lambda*d1x; pts[i1*2+1] += lambda*d1y;
      pts[i2*2+0] += lambda*d2x; pts[i2*2+1] += lambda*d2y;
    }
  }
};

const solve_length_constraints_of_pt_pairs_3d = (pts, i, j, desired_dist, damping) => {
  const [[ x1, y1, z1 ], [ x2, y2, z2 ]] = [ pts.subarray(i*3, i*3+3), pts.subarray(j*3, j*3+3) ];
  const dist = dist3([ x1, y1, z1 ], [ x2, y2, z2 ]);
  const [ d1x, d1y, d1z ] = [ (x1-x2)/dist, (y1-y2)/dist, (z1-z2)/dist ];
  const [ d2x, d2y, d2z ] = [ -d1x, -d1y, -d1z ];
  const grad_sum_sq = d1x**2 + d1y**2 + d1z**2 + d2x**2 + d2y**2 + d2z**2;
  const constraint = dist - desired_dist;
  const lambda = -constraint*damping/grad_sum_sq;
  pts[i*3+0] += lambda*d1x; pts[i*3+1] += lambda*d1y; pts[i*3+2] += lambda*d1z;
  pts[j*3+0] += lambda*d2x; pts[j*3+1] += lambda*d2y; pts[j*3+2] += lambda*d2z;
};

export const solve_length_constraints_tris_3d = ({ pts, tris, orig }, damping=.1) => {
  for (let i = 0; i < tris.length/3; i++) {
    const [ curr_tri, orig_tri, tri_idxs ] = get_tri_with_orig({ pts, orig, tris, dim: 3 }, i);
    const idxs = [ [ 0, 1 ], [ 1, 2 ], [ 2, 0 ] ];
    for (let j = 0; j < idxs.length; j++) {
      const [ j1, j2 ] = idxs[j];
      const [ i1, i2 ] = [ tri_idxs[j1], tri_idxs[j2] ];
      const rest_dist = dist3(orig_tri[j1], orig_tri[j2]);
      solve_length_constraints_of_pt_pairs_3d(pts, i1, i2, rest_dist, damping);
    }
  }
};

export const solve_length_constraints_tris_3d_with_unifrom_dist = ({ pts, tris, rest_dist }, damping=.1) => {
  for (let i = 0; i < tris.length/3; i++) {
    const [ curr_tri, tri_idxs ] = get_tri({ pts, tris, dim: 3 }, i);
    const idxs = [ [ 0, 1 ], [ 1, 2 ], [ 2, 0 ] ];
    for (let j = 0; j < idxs.length; j++) {
      const [ j1, j2 ] = idxs[j];
      const [ i1, i2 ] = [ tri_idxs[j1], tri_idxs[j2] ];
      solve_length_constraints_of_pt_pairs_3d(pts, i1, i2, rest_dist, damping);
    }
  }
};

export const solve_surface_area_constraints_tris_3d = ({ pts, tris, orig }, damping=.1, multiplier=1) => {
  for (let i = 0; i < tris.length/3; i++) {
    const [ curr_tri, orig_tri, tri_idxs ] = get_tri_with_orig({ pts, tris, orig, dim: 3 }, i);
    const tri_area_curr = tri_area(curr_tri);
    const tri_area_orig = tri_area(orig_tri);
    const grads = tri_area_grads(curr_tri);
    let grad_sum_sq = 0;
    for (let j = 0; j < 3; j++) grad_sum_sq += grads[j][0]**2 + grads[j][1]**2 + grads[j][2]**2;
    const constraint = tri_area_curr - tri_area_orig*multiplier;
    const lambda = -constraint*damping/grad_sum_sq;
    for (let j = 0; j < 3; j++) {
      const ii = tri_idxs[j];
      pts[ii*3+0] += lambda*grads[j][0];
      pts[ii*3+1] += lambda*grads[j][1];
      pts[ii*3+2] += lambda*grads[j][2];
    }
  }
};

export const compute_volume_via_stokes_theorem = ({ pts, tris }) => {
  let volume = 0;
  for (let i = 0; i < tris.length/3; i++) {
    const [ [ p1, p2, p3 ], [ i1, i2, i3 ] ] = get_tri({ pts, tris, dim: 3 }, i);
    volume += 1/6 * dot3(p1, cross(p2, p3))
  }

  return volume;
};

export const solve_volume_constraints_tris_3d = ({ pts, tris }, desired_volume, damping=.1, pressure=1) => {
  const grads = Array(pts.length);
  for (let i = 0; i < grads.length; i++) grads[i] = [ 0, 0, 0 ];

  let volume = 0;
  for (let i = 0; i < tris.length/3; i++) {
    const [ [ p1, p2, p3 ], [ i1, i2, i3 ] ] = get_tri({ pts, tris, dim: 3 }, i);
    volume += 1/6 * dot3(p1, cross(p2, p3))
    grads[i1] = add3(grads[i1], cross(p2, p3));
    grads[i2] = add3(grads[i2], cross(p3, p1));
    grads[i3] = add3(grads[i3], cross(p1, p2));
  }

  let grad_sum_sq = 0;
  for (let i = 0; i < grads.length; i++) grad_sum_sq += grads[i][0]**2 + grads[i][1]**2 + grads[i][2]**2;
  const constraint = volume - desired_volume*pressure;
  const lambda = -constraint*damping/grad_sum_sq;

  for (let i = 0; i < tris.length/3; i++) {
    const [ [ p1, p2, p3 ], [ i1, i2, i3 ] ] = get_tri({ pts, tris, dim: 3 }, i);
    pts.set(add3(p1, scale3(lambda, grads[i1])), i1*3);
    pts.set(add3(p2, scale3(lambda, grads[i2])), i2*3);
    pts.set(add3(p3, scale3(lambda, grads[i3])), i3*3);
  }
};

export const apply_distance_based_breakage = ({ pts, tris, orig }, limit_mul=1.5) => {
};

export const make_pbd_boilerplate_from_pts = ({ pts, tris, dim=2 }) => {
  return {
    vels: new Float32Array(pts.length).fill(0),
    prev: pts.slice(), orig: pts.slice(),
    is_fixed: Array(pts.length/dim).fill(false)
  };
};

export const pbd_update = ({ pts, prev, orig, vels, is_fixed, dim=2 },
    { grav=0, num_substeps=1, clamp_box=false, dt=.1 }, constraint_fn=()=>{}, breakage_fn=()=>{}) => {
  const sub_frame_dt = dt/num_substeps;
  for (let itr = 0; itr < num_substeps; itr++) {
    prev.set(pts);
    for (let i = 0; i < pts.length; i++) pts[i] += sub_frame_dt*vels[i];

    constraint_fn();

    if (clamp_box) for (let i = 0; i < pts.length/3; i++) pts[i*3+0] = clamp(pts[i*3+0], -.98, .98);
    if (clamp_box) for (let i = 0; i < pts.length/3; i++) pts[i*3+1] = clamp(pts[i*3+1], -.98, .98);
    for (let i = 0; i < pts.length/dim; i++) if (is_fixed[i]) pts.set(orig.subarray(i*dim, i*dim+dim), i*dim);
    for (let i = 0; i < pts.length; i++) vels[i] = (pts[i] - prev[i])/sub_frame_dt;
    for (let i = 0; i < pts.length/dim; i++) vels[i*dim+1] += grav;
  }

  breakage_fn();
};


export const make_prog = (gl, vsource, fsource) => {
  const prog = gl.createProgram();

  const vert_shader = gl.createShader(gl.VERTEX_SHADER);
  gl.shaderSource(vert_shader, vsource);
  gl.compileShader(vert_shader);
  gl.attachShader(prog, vert_shader);

  const frag_shader = gl.createShader(gl.FRAGMENT_SHADER);
  gl.shaderSource(frag_shader, fsource);
  gl.compileShader(frag_shader);
  gl.attachShader(prog, frag_shader);
  gl.linkProgram(prog);
  gl.useProgram(prog);

  return prog;
};

export const make_gl_index_buf = (gl, arr) => {
  return { arr, type: gl.ELEMENT_ARRAY_BUFFER, buf: gl.createBuffer() };
};

export const make_gl_vert_buf = (gl, arr, loc=0, dim=3) => {
  return { type: gl.ARRAY_BUFFER, buf: gl.createBuffer(), arr, loc, dim };
};

export const draw_indexed_tris = (gl, prog, bufs, [rot1, rot2]=[0, 0]) => {
  let num_tris = 0;
  for (let i = 0; i < bufs.length; i++) {
    if (bufs[i].type == gl.ELEMENT_ARRAY_BUFFER) {
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, bufs[i].buf);
      gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, bufs[i].arr, gl.STATIC_DRAW);
      num_tris = bufs[i].arr.length;
    }

    if (bufs[i].type == gl.ARRAY_BUFFER) {
      gl.bindBuffer(gl.ARRAY_BUFFER, bufs[i].buf);
      gl.enableVertexAttribArray(bufs[i].loc);
      gl.vertexAttribPointer(bufs[i].loc, bufs[i].dim, gl.FLOAT, false, 0, 0);
      gl.bufferData(gl.ARRAY_BUFFER, bufs[i].arr, gl.STREAM_DRAW);
    }
  }

  gl.uniform2f(gl.getUniformLocation(prog, 'rot_ang'), rot1, rot2);
  gl.drawElements(gl.TRIANGLES, num_tris, gl.UNSIGNED_SHORT, 0);
};


//// shaders ////

const diffuse_shading_glsl_fn = () => `
float diffuse_shading(vec3 light, vec3 normal, vec3 pos) {
  vec3 L = light - pos;
  //float NdotL = max(dot(L, normal), dot(L, -normal));
  float ambient = .3;
  float NdotL = dot(L, normal);
  float diffuse = max(0., NdotL);
  return diffuse*.3 + ambient;
}`;

const rotation_glsl_fn = () => `
  uniform vec2 rot_ang;
  vec3 rotate(vec3 pos) {
    mat3 rot_xz = mat3(
      cos(rot_ang.x), 0., -sin(rot_ang.x),
      0., 1., 0.,
      sin(rot_ang.x), 0.,  cos(rot_ang.x)
    );
    mat3 rot_yz = mat3(
      1., 0., 0.,
      0., cos(rot_ang.y), -sin(rot_ang.y),
      0., sin(rot_ang.y),  cos(rot_ang.y)
    );
    return rot_yz * rot_xz * pos;
  }
`;

const rnd_glsl_fn = () => `
  float rand(vec2 co){ return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453); }`;

// only works for some triangle setups (for some reason)
export const wireframe_shader = (gl, line_width=10) => [ `#version 300 es
  layout (location=0) in vec3 pos;
  out vec3 vpos;
  ${rotation_glsl_fn()}
  void main() {
    vec3 vert_pos = pos*vec3(1, -1, 1);
    gl_Position = vec4(rotate(vert_pos), 1);
    if (gl_VertexID % 3 == 0) vpos = vec3(1, 0, 0);
    if (gl_VertexID % 3 == 1) vpos = vec3(0, 1, 0);
    if (gl_VertexID % 3 == 2) vpos = vec3(0, 0, 1);
  }`, `#version 300 es
  precision highp float;
  in vec3 vpos;
  out vec4 frag_color;
  void main() {
    frag_color = vec4(0, 0, 0, 0);
    float th = 1./float(${gl.canvas.width}) * float(${line_width});
    if (vpos.x < th || vpos.y < th || vpos.z < th) frag_color = vec4(.9, .9, .9, 1);
  }`
];

export const normals_shader = () => [ `#version 300 es
  layout (location=0) in vec3 pos;
  layout (location=1) in vec3 normal;
  out vec3 vnormal;
  ${rotation_glsl_fn()}
  void main() {
    vec3 p = rotate(pos*vec3(1, -1, 1));
    gl_Position = vec4(p, 1); vnormal = normal;
  }`, `#version 300 es
  precision highp float;
  in vec3 vnormal;
  out vec4 frag_color;
  void main() { frag_color = vec4(abs(vnormal), 1);
  //void main() { frag_color = vec4(max(0., vnormal.z), max(0., -vnormal.z), 0, 1);
  }`
];

export const diffuse_lighting_shader = () => [ `#version 300 es
  layout (location=0) in vec3 pos;
  layout (location=1) in vec3 normal;
  out float vshading;
  ${rotation_glsl_fn()}
  ${diffuse_shading_glsl_fn()}
  void main() {
    vec3 vert_pos = pos*vec3(1, -1, 1);
    vec3 vpos = rotate(vert_pos.xyz);
    gl_Position = vec4(vpos, 1);

    vec3 light = vec3(0, -2, -2);
    vshading = diffuse_shading(light, normal, vert_pos);
  }`, `#version 300 es
  precision highp float;
  in float vshading;
  out vec4 frag_color;
  void main() {
    vec3 color = vec3(.2, .2, .9);
    frag_color = vec4(color*vshading, 1);
  }`
];

export const curvature_shader = () => [ `#version 300 es
  layout (location=0) in vec3 pos;
  layout (location=1) in vec3 normal;
  //out float vshading;
  out vec3 vnormal;
  out vec3 vpos;
  ${rotation_glsl_fn()}
  void main() {
    vec3 vert_pos = pos*vec3(1, -1, -1);
    vpos = rotate(vert_pos.xyz);
    gl_Position = vec4(vpos, 1);
    vnormal = normal;
  }`, `#version 300 es
  precision highp float;
  //in float vshading;
  in vec3 vpos;
  in vec3 vnormal;
  out vec4 frag_color;
  void main() { // via evanw
    vec3 n = normalize(vnormal);
    vec3 dx = dFdx(n);
    vec3 dy = dFdy(n);
    vec3 xneg = n - dx;
    vec3 xpos = n + dx;
    vec3 yneg = n - dy;
    vec3 ypos = n + dy;
    float depth = length(vpos);
    float curvature = (cross(xneg, xpos).y - cross(yneg, ypos).x) * 4.0 / depth;
    frag_color = vec4(vec3(curvature + .5), 1);
  }`
];


export const rnd_color_shader = () => [ `#version 300 es
  layout (location=0) in vec3 pos;
  flat out vec3 vcol;
  out vec3 vpos;
  ${rnd_glsl_fn()}
  ${rotation_glsl_fn()}
  void main() {
    vec3 p = rotate(pos*vec3(1, -1, 1));
    gl_Position = vec4(p, 1);
    vpos = p;
    vcol = vec3(
      rand(vec2(gl_VertexID, gl_VertexID*2)),
      rand(vec2(gl_VertexID, gl_VertexID*3)),
      rand(vec2(gl_VertexID, gl_VertexID*4))
    );
  }`, `#version 300 es
  precision highp float;
  flat in vec3 vcol;
  in vec3 vpos;
  out vec4 frag_color;
  void main() {
    float shading = 1./(vpos.z+1.2);
    frag_color = vec4(vcol * shading, 1);
  }`
];

export const screenspace_lighting_shader = () => [ `#version 300 es
  layout (location=0) in vec3 pos;
  out vec3 vpos;
  ${rotation_glsl_fn()}
  void main() {
    vec3 vert_pos = pos*vec3(1, -1, 1);
    gl_Position = vec4(rotate(vert_pos.xyz), 1);
    vpos = vert_pos;
  }`, `#version 300 es
  precision highp float;
  in vec3 vpos;
  out vec4 frag_color;
  ${diffuse_shading_glsl_fn()}
  void main() {
    vec3 light = vec3(0, -2, 2);
    vec3 color = vec3(.2, .2, .9);
    vec3 normal = normalize(cross(dFdx(vpos), dFdy(vpos)));
    frag_color = vec4(color*diffuse_shading(light, normal, vpos), 1);
  }`
];


