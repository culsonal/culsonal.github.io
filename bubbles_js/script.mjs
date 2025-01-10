
import * as lib from './lib.mjs';
import * as meshes from './meshes.mjs';

const gl  = document.getElementById('c').getContext('webgl2'); 
gl.enable(gl.DEPTH_TEST);

let [ pts, tris ] = meshes.sphere2;
[ pts, tris ] = lib.subdivide_and_normalize({ pts, tris });
[ pts, tris ] = lib.subdivide_and_normalize({ pts, tris }, .4);
[ pts, tris ] = lib.subdivide_and_normalize({ pts, tris }, .4);

const geom = {
  pts, tris, normals: lib.make_normals({ pts, tris }), dim: 3,
  ...lib.make_pbd_boilerplate_from_pts({ pts, tris, dim: 3 }),
};

const gl_bufs = [
  lib.make_gl_index_buf(gl, tris),
  lib.make_gl_vert_buf(gl, pts, 0),
];

const desired_volume = lib.compute_volume_via_stokes_theorem({ pts, tris });

const [ vsrc, fsrc ] = lib.screenspace_lighting_shader();
const prog = lib.make_prog(gl, vsrc, fsrc);

let itr = 0;
const loop = () => {
  window.requestAnimationFrame(loop);
  itr += 1;
  lib.draw_indexed_tris(gl, prog, gl_bufs, [ 0, 3.1415/4 ]);

  lib.pbd_update(geom, { grav: .001, clamp_box: true, num_substeps: 1 }, () => {
    lib.solve_length_constraints_tris_3d(geom, .01);
    lib.solve_volume_constraints_tris_3d(geom, desired_volume, .006, 1.1);
    lib.solve_surface_area_constraints_tris_3d(geom, 1, .7);
  });
};

loop();

