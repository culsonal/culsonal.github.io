
import * as gpu from './gpu.mjs';

const gl = document.getElementById('c').getContext('webgl2');
gl.getExtension("EXT_color_buffer_float");
gl.canvas.style.width = gl.canvas.style.height = 800;
gl.canvas.width = gl.canvas.height = 3000;
gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

const draw_pts = gpu.make_prog(gl, `${await (await fetch('hash.c')).text()}
  uniform vec2 ellipse_params;
  uniform vec2 obj_center;
  uniform float ior;
  uniform int itr;

  vec2 sample_half_circle(vec2 normal, vec2 seed) {
    int attempt = 0;
    vec2 v;
    do {
      v = hash23(vec3(seed, attempt))*2.-1.;
      attempt += 1;
    } while (dot(v, normal) < 0.);
    return normalize(v);
  }

  float implicit_ellipse(vec2 pt, vec2 ab) {
    return pow(pt.x/ab.x, 2.) + pow(pt.y/ab.y, 2.) - 1.;
  }

  float implicit_heart(vec2 pt) {
    return pow(pt.x*pt.x+pt.y*pt.y-1., 3.) - pt.x*pt.x * pt.y*pt.y*pt.y;
  }

  float implicit_union(float f1, float f2) {
    return min(f1, f2);
  }

  float implicit_fn(vec2 pt) {
    //return implicit_ellipse(pt-obj_center, ellipse_params);
    return implicit_union(
      implicit_heart((pt-obj_center)*3.5),
      implicit_ellipse(pt, vec2(.2, .2))
    );
  }

  bool is_inside_impl_fn(vec2 pt) {
    return implicit_fn(pt) < 0.;
  }

  vec2 implicit_normal(vec2 pt) {
    const float eps = .0001;
    return normalize(vec2(
      (implicit_fn(pt+vec2(eps, 0)) - implicit_fn(pt-vec2(eps, 0)))/(2.*eps),
      (implicit_fn(pt+vec2(0, eps)) - implicit_fn(pt-vec2(0, eps)))/(2.*eps)
    ));
  }

  void main() {
    vec2 seed = vec2(itr, gl_VertexID);
    vec2 light_a = vec2(-.2, .9);
    vec2 light_b = vec2( .2, .9);
    vec2 light_normal = (light_a-light_b).yx * vec2(-1, 1);
    float m = hash13(vec3(seed, 1));
    vec2 pos = mix(light_a, light_b, m);
    vec2 dir = sample_half_circle(light_normal, seed);

    const float max_steps = 3000.;
    int num_steps = int(max_steps * hash13(vec3(seed, 2)));
    const float step_size = .0008;
    for (int i = 0; i < num_steps; i++) {
      vec2 pos_new = pos + dir*step_size;

      if (is_inside_impl_fn(pos_new) && !is_inside_impl_fn(pos)) {
	vec2 en = implicit_normal(pos_new);
	dir = normalize(refract(vec3(dir, 0), vec3(en, 0), ior).xy);
      }
      if (!is_inside_impl_fn(pos_new) && is_inside_impl_fn(pos)) {
	vec2 en = implicit_normal(pos_new);
	dir = normalize(refract(vec3(dir, 0), vec3(en, 0), 1./ior).xy);
      }

      if (pos_new.x >  .95) dir.x *= -1.;
      if (pos_new.x < -.95) dir.x *= -1.;

      pos = pos_new;
    }

    gl_Position = vec4(pos, 0, 1);
    gl_PointSize = 4.;
  }`, `
  out vec4 col;
  void main() { col = vec4(vec3(.7), 1); }`);

let [ mx, my ] = [ .5, .5 ];
gl.canvas.onmousemove = e => [ mx, my ] = [ e.offsetX/800, 1-e.offsetY/800 ];
let itr = 0; const loop = () => { window.requestAnimationFrame(loop); itr += 1;

gl.useProgram(draw_pts);
gl.uniform1i(gl.getUniformLocation(draw_pts, 'itr'), itr);
gl.uniform2f(gl.getUniformLocation(draw_pts, 'ellipse_params'),
  Math.sin(itr/50)*.35+.4,
  Math.cos(itr/100+.2)*.3+.4);
gl.uniform1f(gl.getUniformLocation(draw_pts, 'ior'), .6);//.4*Math.sin(itr/50)+.5);
gl.uniform2f(gl.getUniformLocation(draw_pts, 'obj_center'), mx*2-1, my*2-1);
gl.drawArrays(gl.POINTS, 0, 100_000);

}; loop();
