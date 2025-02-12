
import * as gpu from './lib/gpu.mjs';

const canvas = document.getElementById('c');
const gl = gpu.good_defaults_gl(canvas, false, true);

const prog = gpu.make_frag_prog(gl, gl.canvas.width, `
  uniform vec3 index_of_refraction;
  uniform vec3 camera_xyz;
  uniform vec3 sphere_xyz;

  vec3 background(vec2 xy) {
    float freq = 10.;
    float x = xy.x*.5+1000.;
    float y = xy.y*.5+1000.;
    if (int(x * freq) % 2 == 0 && int(y * freq) % 2 == 0) return vec3(1);
    if (int(x * freq) % 2 == 1 && int(y * freq) % 2 == 1) return vec3(1);
    return vec3(0);
  }

  bool is_inside_sphere(vec3 pt) {
    return length(pt - sphere_xyz) < 1.2;
  }

  vec3 sphere_normal(vec3 pt) {
    return normalize(pt - sphere_xyz);
  }

  vec3 get_color(float ior) {
    vec2 xy = gl_PointCoord.xy*2.-1.;
    vec3 pt_on_camera_plane = vec3(xy, -2) + (camera_xyz+vec3(0, 0, 3));
    vec3 camera_pt = camera_xyz; // vec3(0, 0, -3);
    vec3 ray_dir = normalize(pt_on_camera_plane - camera_pt);
    float eps_move = .01;
    int max_num_iters = 1000;
    vec3 pt = camera_pt;
    vec3 color = vec3(0);
    for (int i = 0; i < max_num_iters; i++) {
      vec3 pt_new = pt+ray_dir*eps_move;
      if (is_inside_sphere(pt) && !is_inside_sphere(pt_new) ||
         !is_inside_sphere(pt) &&  is_inside_sphere(pt_new)) {
        ray_dir = normalize(refract(ray_dir, sphere_normal(pt), ior));
      }

      pt = pt_new;

      if (pt_new.z > 1.1) {
        color = background(pt_new.xy);
        break;
      }
    }

    return color;
  }

  out vec3 frag_color;
  void main() {
    float r = get_color(index_of_refraction.r).r;
    float g = get_color(index_of_refraction.g).g;
    float b = get_color(index_of_refraction.b).b;

    frag_color = vec3(r, g, b);
    //frag_color = background(xy);
  }`);

let itr = 0; const loop = () => { itr += 1; window.requestAnimationFrame(loop);

  //const index_of_refraction = (Math.sin(itr/100)*.5+.5)*.6 + .05;
  //const index_of_refraction = .35;

  gpu.run_prog_pts(gl, prog, 1, {
    //index_of_refraction: [ .3+Math.sin(itr/15)*.008, .3+Math.cos(itr/20)*.005, .3+Math.sin(itr/10)*.005 ],
    index_of_refraction: [ .3+.008, .3-.005, .3 ],
    //camera_xyz: [ 0, 0, -3 ],
    camera_xyz: [ Math.sin(itr/200)*.4, Math.cos(itr/150)*.4, -3 ],
    sphere_xyz: [ Math.sin(itr/20)*.2, Math.cos(itr/15)*.2, 0 ]
  });
  //console.log(index_of_refraction);

}; loop();

