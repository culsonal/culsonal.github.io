
// wasm lib

import wasm_ from './lib_wasm.mjs';
export const wasm = await wasm_();

const wasm_mem = wasm.asm.b.buffer;

export const new_float32array = (size=-1, data=null) => {
  const ptr = wasm._alloc_float32array(size);
  const arr = new Float32Array(wasm_mem, ptr, size);
  if (data != null) arr.set(data);
  return { ptr, size, arr };
};

export const new_uint8array = (size=-1, data=null) => {
  const ptr = wasm._alloc_uint8array(size);
  const arr = new Uint8ClampedArray(wasm_mem, ptr, size);
  if (data != null) arr.set(data);
  return { ptr, size, arr };
};

export const new_uint16array = (size=-1, data=null) => {
  const ptr = wasm._alloc_uint16array(size);
  const arr = new Uint16Array(wasm_mem, ptr, size);
  if (data != null) arr.set(data);
  return { ptr, size, arr };
};

export const new_uint32array = (size=-1, data=null) => {
  const ptr = wasm._alloc_uint32array(size);
  const arr = new Uint32Array(wasm_mem, ptr, size);
  if (data != null) arr.set(data);
  return { ptr, size, arr };
};

export const fill_arr_pt_cloud_3d_sphere = (pts, radius=1) =>
  wasm._fill_arr_pt_cloud_3d_sphere(pts.ptr, pts.size/3, radius);
export const fill_arr_pt_cloud_2d_cube = (pts, cube_size=1) =>
  wasm._fill_arr_pt_cloud_2d_cube(pts.ptr, pts.size/2, cube_size);
export const fill_arr_pt_cloud_3d_cube = (pts, cube_size=1) =>
  wasm._fill_arr_pt_cloud_3d_cube(pts.ptr, pts.size/3, cube_size);
export const move_pts_2d = ({ pts, num_pts=1 }, [ox, oy]=[0, 0]) => wasm._move_pts_2d(pts.ptr, num_pts, ox, oy);
export const move_pts_3d = ({ pts, num_pts=1 }, [ox, oy, oz]=[0, 0, 0]) =>
  wasm._move_pts_3d(pts.ptr, num_pts, ox, oy, oz);
export const move_pts = ({ pts, num_pts=1, dim=2 }, offset) => dim == 2 ?
  move_pts_2d(pts.ptr, num_pts, offset[0], offset[1]) : move_pts_3d(pts.ptr, num_pts, offset[0], offset[1], offset[2]);

export const modify_vels_via_curl_noise_2d = (a, num_basis_fns=5, strength=1, noise_seed=3) =>
  wasm._modify_vels_via_curl_noise_2d(a.pts.ptr, a.num_pts, a.inv_mass, a.vels.ptr, a.dt, num_basis_fns,
    strength, noise_seed);

export const pbd_pre = ({ pts, prev, vels, num_pts=1, dt=.1, num_substeps=1, dim=2 }) =>
  wasm._pbd_pre(pts.ptr, num_pts, prev.ptr, vels.ptr, dt, num_substeps, dim);
export const pbd_post = ({ pts, prev, vels, orig, inv_mass, num_pts=1, dt=.1, num_substeps=1, dim=2, grav=0,
  particle_radius=.03, air_mass=0, vel_damping=1, clamp_val=.98, clamp_box=true, limit_max_vel=true }) =>
    wasm._pbd_post(pts.ptr, num_pts, prev.ptr, vels.ptr, orig.ptr, inv_mass.ptr, dt, num_substeps, dim, grav,
      particle_radius, air_mass, vel_damping, clamp_val, clamp_box, limit_max_vel);

export const make_pbd_boilerplate_from_pts = ({ pts, num_pts, dim }) => {
  const inv_mass = new_float32array(num_pts*dim, pts.arr);
  inv_mass.arr.fill(1);
  const is_pt_active = new_uint8array(num_pts);
  is_pt_active.arr.fill(1);
  return {
    prev: new_float32array(num_pts*dim, pts.arr),
    orig: new_float32array(num_pts*dim, pts.arr),
    vels: new_float32array(num_pts*dim),
    is_pt_active, inv_mass
  };
};

export const make_pbd_boilerplate_from_pts_and_tris = ({ pts, num_pts, num_tris, dim }) => {
  const is_tri_active = new_uint8array(num_tris);
  is_tri_active.arr.fill(1);
  return { is_tri_active, ...make_pbd_boilerplate_from_pts({ pts, num_pts, dim }) };
};

export const make_img_boilerplate = ctx => {
  const { ptr, arr, size } = new_uint8array(ctx.canvas.width*ctx.canvas.height*4);
  return {
    data: new ImageData(arr, ctx.canvas.width, ctx.canvas.height),
    depth: lib.new_float32array(ctx.canvas.width*ctx.canvas.height),
    ptr, arr, size, img_width: ctx.canvas.width, img_height: ctx.canvas.height
  };
};

export const make_hash_boilerplate = num_pts => {
  return { hash: new_uint16array(num_pts*3), pts_reordered: new_uint16array(num_pts) } };

export const compute_hash = a =>
  wasm._compute_hash(a.pts.ptr, a.num_pts, a.pts_reordered.ptr, a.hash.ptr, a.particle_radius, a.dim);

export const draw_surface_pts_2d_cpu = (ctx, a, img) =>
  { wasm._draw_surface_pts_2d(a.pts.ptr, a.num_pts, a.pts_reordered.ptr, a.hash.ptr, a.particle_radius,
    img.ptr, img.img_width, img.img_height, a.pt_size); ctx.putImageData(img.data, 0, 0); };

export const draw_pts_2d_cpu = (ctx, a, img) =>
  { wasm._draw_pts_2d(a.pts.ptr, a.num_pts, img.ptr, img.img_width, img.img_height, a.pt_size);
    ctx.putImageData(img.data, 0, 0); };

export const draw_pts_3d_cpu = (ctx, a, img) => {
  wasm._draw_pts_3d(a.pts.ptr, a.num_pts, img.ptr, img.depth.ptr, img.img_width, img.img_height, a.pt_size,
    a.offset[0], a.offset[1], a.offset[2], a.rot[0], a.rot[1], a.rot[2], a.max_depth);
  ctx.putImageData(img.data, 0, 0); };

export const draw_pts_cpu = (ctx, a, img) => a.dim == 2 ?
  draw_pts_2d_cpu(ctx, a, img) : draw_pts_3d_cpu(ctx, a, img);

export const solve_pairwise_constraints = (a, damping=.9, sticky_mul=1) => a.dim == 2 ?
  wasm._solve_pairwise_constraints_2d(a.pts.ptr, a.num_pts, a.pts_reordered.ptr, a.inv_mass.ptr, a.hash.ptr,
    a.particle_radius, damping, sticky_mul) :
  wasm._solve_pairwise_constraints_3d(a.pts.ptr, a.num_pts, a.pts_reordered.ptr, a.inv_mass.ptr, a.hash.ptr,
    a.particle_radius, damping, sticky_mul);

export const solve_neohookean_constraints_tris_2d = (a, [d1, d2]) => wasm._solve_neohookean_constraints_tris_2d(
  a.pts.ptr, a.num_pts, a.tris.ptr, a.num_tris, a.orig.ptr, a.is_tri_active.ptr, a.inv_mass.ptr, d1, d2);

export const solve_length_constraints_tris_2d = (a, damping) => wasm._solve_length_constraints_tris_2d(a.pts.ptr,
  a.num_pts, a.tris.ptr, a.num_tris, a.orig.ptr, a.is_tri_active.ptr, a.inv_mass.ptr, damping);

export const solve_length_constraints_tris_3d = (a, damping) => wasm._solve_length_constraints_tris_3d(a.pts.ptr,
  a.num_pts, a.tris.ptr, a.num_tris, a.orig.ptr, a.is_tri_active.ptr, a.inv_mass.ptr, damping);

export const solve_area_constraints_tris_2d = (a, damping) => wasm._solve_area_constraints_tris_2d(a.pts.ptr,
  a.num_pts, a.tris.ptr, a.num_tris, a.orig.ptr, a.is_tri_active.ptr, a.inv_mass.ptr, damping);

export const solve_surface_area_constraints_tris_3d = (a, dmp) => wasm._solve_surface_area_constraints_tris_2d(
  a.pts.ptr, a.num_pts, a.tris.ptr, a.num_tris, a.orig.ptr, a.is_tri_active.ptr, a.inv_mass.ptr, dmp);

export const solve_volume_constraints_tris_3d = (a, dmp=.1, pressure=1) =>
  wasm._solve_volume_constraints_tris_3d(a.pts.ptr, a.num_pts, a.tris.ptr, a.num_tris, a.orig.ptr, a.is_tri_active.ptr,
    a.inv_mass.ptr, dmp, pressure);

export const apply_energy_based_breakage_tris_2d = (a, [ thr1, thr2 ]) =>
  wasm._apply_energy_based_breakage_tris_2d(a.pts.ptr, a.num_pts, a.tris.ptr, a.num_tris, a.orig.ptr,
    a.is_tri_active.ptr, thr1, thr2);

export const apply_distance_based_breakage_tris_2d = (a, threshold) =>
  wasm._apply_distance_based_breakage_tris_2d(a.pts.ptr, a.num_pts, a.tris.ptr, a.num_tris, a.orig.ptr,
    a.is_tri_active.ptr, threshold);

export const apply_distance_based_breakage_tris_3d = (a, threshold) =>
  wasm._apply_distance_based_breakage_tris_3d(a.pts.ptr, a.num_pts, a.tris.ptr, a.num_tris, a.orig.ptr,
    a.is_tri_active.ptr, threshold);

// js lib

export const sleep = ms => new Promise(r => setTimeout(r, ms));
export const dist2 = ([ x1, y1 ], [ x2, y2 ]) => Math.sqrt((x1-x2)**2 + (y1-y2)**2);
export const dist3 = ([ x1, y1, z1 ], [ x2, y2, z2 ]) => Math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2);
export const dist  = (p1, p2) => p1.length == 2 ? dist2(p1, p2) : dist3(p1, p2);
export const len2    = ([ x, y ]) => Math.sqrt(x*x + y*y);
export const len2_sq = ([ x, y ]) => x*x + y*y;
export const len3    = ([ x, y, z ]) => Math.sqrt(x*x + y*y + z*z);
export const len3_sq = ([ x, y, z ]) => x*x + y*y + z*z;
export const len    = pt => pt.length == 2 ? len2(pt) : len3(pt);
export const len_sq = pt => pt.length == 2 ? len_sq2(pt) : len_sq3(pt);
export const normalize2 = ([ x, y ]) => { const l = len2([ x, y ]); return [ x/l, y/l ]; };
export const normalize3 = ([ x, y, z ]) => { const l = len3([ x, y, z ]); return [ x/l, y/l, z/l ]; };
export const normalize = pt => pt.length == 2 ? normalize2(pt) : normalize3(pt);
export const add2 = ([ x1, y1 ], [ x2, y2 ]) => [ x1+x2, y1+y2 ];
export const add3 = ([ x1, y1, z1 ], [ x2, y2, z2 ]) => [ x1+x2, y1+y2, z1+z2 ];
export const add  = (p1, p2) => p1.length == 2 ? add2(p1, p2) : add3(p1, p2);
export const sub2 = ([ x1, y1 ], [ x2, y2 ]) => [ x1-x2, y1-y2 ];
export const sub3 = ([ x1, y1, z1 ], [ x2, y2, z2 ]) => [ x1-x2, y1-y2, z1-z2 ];
export const sub  = (p1, p2) => p1.length == 2 ? sub2(p1, p2) : sub3(p1, p2);
export const dot2 = ([ x1, y1 ], [ x2, y2 ]) => x1*x2 + y1*y2;
export const dot3 = ([ x1, y1, z1 ], [ x2, y2, z2 ]) => x1*x2 + y1*y2 + z1*z2;
export const dot  = (p1, p2) => p1.length == 2 ? dot2(p1, p2) : dot3(p1, p2);
export const mix2 = ([ x1, y1 ], [ x2, y2 ], m) => [ x1*m+(1-m)*x2, y1*m+(1-m)*y2 ];
export const mix3 = ([ x1, y1, z1 ], [ x2, y2, z2 ], m) => [ x1*m+(1-m)*x2, y1*m+(1-m)*y2, z1*m+(1-m)*z2 ];
export const mix = (p1, p2, m) => p1.length == 2 ? mix2(p1, p2, m) : mix3(p1, p2, m);
export const cross = ([ v0, v1, v2 ], [ w0, w1, w2 ]) => [ v1*w2-v2*w1, v2*w0-v0*w2, v0*w1-v1*w0 ];
export const clamp = (v, min, max) => Math.max(Math.min(v, max), min);
export const clamp2 = ([ x, y ], min, max) => [ clamp(x, min, max), clamp(y, min, max) ];
export const orth = ([ x, y ]) => [ -y, x ];
export const eye_matr3 = () => [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
export const mat_vec_mul2 = ([m1, m2], v) => [dot2(m1, v), dot2(m2, v)];
export const mat_vec_mul3 = ([m1, m2, m3], v) => [dot3(m1, v), dot3(m2, v), dot3(m3, v)];
export const scale2 = (s, [x, y]) => [s*x, s*y];
export const scale3 = (s, [x, y, z]) => [s*x, s*y, s*z];
export const scale  = (s, pt) => pt.length == 2 ? scale2(s, pt) : scale3(s, pt);
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
export const avg = pts => {
  let avg = Array(pts[0].length).fill(0);
  for (let i = 0; i < pts.length; i++) avg = add(avg, pts[i]);
  return scale(1/pts.length, avg);
};

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

export const tri_area_2d_and_grads = ([[a, b], [c, d], [e, f]]) => {
  const signed_area = (c-a)*(f-b) - (d-b)*(e-a);
  const coeff = 1/2 * Math.sign(signed_area);
  const dp1 = [ coeff*(d-f), coeff*(e-c) ];
  const dp2 = [ coeff*(f-b), coeff*(a-e) ];
  const dp3 = [ coeff*(b-d), coeff*(c-a) ];
  return [ Math.abs(signed_area), [ dp1, dp2, dp3 ] ];
};

// via sympy
const tri_area_grads_3d = ([ p1, p2, p3 ]) => {
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

export const good_defaults_ctx = canvas => canvas.getContext('2d', {
  alpha: false, desynchronized: true, antialias: false,
  powerPreference: 'high-performance', preserveDrawingBuffer: true,
  colorSpace: 'display-p3'
});

export const good_defaults_gl = canvas => {
  const gl = canvas.getContext('webgl2'); 
  gl.enable(gl.DEPTH_TEST);
  //gl.enable(gl.CULL_FACE);
  //gl.cullFace(gl.FRONT_AND_BACK); 
  //gl.cullFace(gl.FRONT); gl.cullFace(gl.BACK); 
  return gl;
};

export const isometric_bending_grads = ([ p1, p2, p3, p4 ]) => {};

export const empty_2d_array = (dim1, dim2) => {
  const arr = Array(dim1);
  for (let i = 0; i < dim1; i++) arr[i] = Array(dim2).fill(0);
  return arr;
};

export const draw_tris_2d = (ctx, pts, tris, is_active, lwidth=1, color='#fff') => {
  const [ w, h ] = [ ctx.canvas.width, ctx.canvas.height ];
  ctx.lineWidth = lwidth; ctx.strokeStyle = color;
  for (let i = 0; i < tris.length/3; i++) {
    if (is_active[i] == 0) continue;
    const [[ x1, y1 ], [ x2, y2 ], [ x3, y3 ]] = get_tri({ pts, tris, dim: 2 }, i)[0];
    ctx.beginPath();
    ctx.moveTo((x1*.5+.5)*w, (y1*.5+.5)*h);
    ctx.lineTo((x2*.5+.5)*w, (y2*.5+.5)*h);
    ctx.lineTo((x3*.5+.5)*w, (y3*.5+.5)*h);
    ctx.closePath();
    ctx.stroke();
  }
};

export const draw_pts_2d = (ctx, pts, size=2, color='#fff') => {
  const [ w, h ] = [ ctx.canvas.width, ctx.canvas.height ];
  ctx.fillStyle = color;
  for (let i = 0; i < pts.length/2; i++) {
    const [ x, y ] = pts.subarray(i*2, i*2+2);
    ctx.fillRect((x*.5+.5)*w-size/2, (y*.5+.5)*h-size/2, size, size);
  }
};

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
};

export const make_normals = ({ pts, tris }) => {
  const normals = new Float32Array(pts.length);
  update_normals({ pts, tris, normals });
  return normals;
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
  return { arr, type: gl.ARRAY_BUFFER, buf: gl.createBuffer(), loc, dim };
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
    return rot_xz * rot_yz * pos;
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
  void main() { frag_color = vec4(max(0., vnormal.z), max(0., -vnormal.z), 0, 1);
  //void main() { frag_color = vec4(abs(vnormal), 1);
  //void main() { frag_color = vec4(max(0., vnormal.z), max(0., -vnormal.z), 0, 1);
  }`
];

export const diffuse_lighting_shader = () => [ `#version 300 es
  layout (location=0) in vec3 pos;
  layout (location=1) in vec3 normal;
  out float vshading;
  out vec3 vpos;
  ${rotation_glsl_fn()}
  ${diffuse_shading_glsl_fn()}
  void main() {
    vec3 vert_pos = pos*vec3(1, -1, 1);
    vpos = rotate(vert_pos.xyz);
    gl_Position = vec4(vpos.xy+vec2(0, .5), 0, vpos.z+.5);

    vec3 light = vec3(0, -2, -2);
    vshading = diffuse_shading(light, normalize(normal), vert_pos);
  }`, `#version 300 es
  precision highp float;
  in float vshading;
  in vec3 vpos;
  out vec4 frag_color;
  void main() {
    vec3 color = vec3(.2, .2, .9);
    frag_color = vec4(color*vshading, 1);
    gl_FragDepth = vpos.z/10.;
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
  precision lowp float;
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


