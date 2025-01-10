
import * as lib from './lib/lib.mjs';
import * as meshes from './meshes.mjs';

const dim = 3;

const gl = lib.good_defaults_gl(document.getElementById('c'));
gl.canvas.width  = Math.min(window.innerWidth, window.innerHeight)/1;
gl.canvas.height = Math.min(window.innerWidth, window.innerHeight)/1;
gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
gl.canvas.style = 'height:100%;width:auto';

const sphere_pts = meshes.sphere10.pts.threejs.slice();
const sphere_tris = meshes.sphere10.tris.slice();

const pts = lib.new_float32array(sphere_pts.length*2);
const tris = lib.new_uint16array(sphere_tris.length*2);
pts.arr.set(sphere_pts, 0); tris.arr.set(sphere_tris, 0);
pts.arr.set(sphere_pts.slice(), sphere_pts.length); tris.arr.set(sphere_tris.slice(), sphere_tris.length);
for (let i = sphere_tris.length; i < sphere_tris.length*2; i++) tris.arr[i] += sphere_pts.length/dim;

const num_pts  = pts.arr.length/dim;
const num_tris = tris.arr.length/3;

for (let i = 0; i < num_pts*dim; i++) pts.arr[i] *= .2;
for (let i = Math.floor(num_pts/2); i < num_pts; i++) pts.arr[i*dim+1] += -.5;

const particle_radius = .031;

const sim = {
  pts, tris, num_pts, num_tris, dim, particle_radius, grav: .001,
  num_substeps: 3, vel_damping: .999, dt: .1,
  pt_size: 2, clamp_val: .9, offset: [0, 0, 1], rot: [0, 0, .2], max_depth: 4., // 3d
  ...lib.make_pbd_boilerplate_from_pts_and_tris({ pts, num_pts, tris, num_tris, dim }),
  ...lib.make_hash_boilerplate(num_pts),
};

console.log('num pts', sim.num_pts);

const gl_bufs = [ lib.make_gl_index_buf(gl, sim.tris.arr), lib.make_gl_vert_buf(gl, sim.pts.arr, 0) ];
const [ vsrc, fsrc ] = lib.screenspace_lighting_shader();
const prog = lib.make_prog(gl, vsrc, fsrc);

let itr = 0; const loop = () => { itr += 1;
  window.requestAnimationFrame(loop);
  lib.draw_indexed_tris(gl, prog, gl_bufs, [ 0, .2 ]);

  lib.compute_hash(sim);
  for (let i = 0; i < sim.num_substeps; i++) {
    lib.pbd_pre(sim);
    lib.solve_pairwise_constraints(sim, 1);
    lib.solve_length_constraints_tris_3d(sim, .35);
    lib.solve_volume_constraints_tris_3d(sim, .2, 1);
    lib.pbd_post(sim);
  }
};

loop();

