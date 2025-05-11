
import * as gpu from './gpu.mjs';

const wasm = (await WebAssembly.instantiateStreaming(fetch('a.out.wasm'))).instance.exports;
const dim = wasm.get_dim();
const tris = new Float32Array(wasm.memory.buffer, wasm.get_tris_ptr(), 3*2*3*dim**2);
const colors = new Float32Array(wasm.memory.buffer, wasm.get_colors_ptr(), 3*2*4*dim**2);
const normals = new Float32Array(wasm.memory.buffer, wasm.get_normals_ptr(), 3*2*3*dim**2);
const gl = document.getElementById('c').getContext('webgl2'); gl.enable(gl.DEPTH_TEST);
gl.canvas.style.width = gl.canvas.style.height = Math.min(innerWidth, innerHeight)-2;
gl.canvas.width = gl.canvas.height = parseInt(gl.canvas.style.width)/2;
//let [ mx, my, mvx, mvy, mouse_active ] = [ .5, .5, 0, 0, 0 ];
let [ mx, my, mvx, mvy, mouse_active ] = [ .52, .18, 0, 0, 0 ];
gl.canvas.onpointerdown = e => mouse_active = 1;
gl.canvas.onpointerup   = e => mouse_active = 0;
gl.canvas.onpointermove = e => [ mvx, mvy, mx, my ] = [
  e.movementX/parseInt(gl.canvas.style.width), 1-e.movementY/parseInt(gl.canvas.style.height),
  e.offsetX/parseInt(gl.canvas.style.width), 1-e.offsetY/parseInt(gl.canvas.style.height) ];

const perlin_scale = 1.5;
const octaves = 5;
const persitence = .6;

const seed = Math.random()*10000;
wasm.init(seed, perlin_scale, octaves, persitence); console.log('seed', Math.floor(seed));

const draw = gpu.make_prog(gl, `
  uniform vec3 mouse;
  uniform float itr;
  layout (location=0) in vec3 pos;
  layout (location=1) in vec4 color;
  layout (location=2) in vec3 normal;
  out float vitr;
  out vec3 vpos;
  out vec4 vcolor;
  out vec3 vnormal;

  vec3 rotate(vec3 pos, vec2 rot) {
    return mat3(cos(rot.x), 0., -sin(rot.x), 0., 1., 0., sin(rot.x), 0., cos(rot.x)) *
	   mat3(1., 0., 0., 0., cos(rot.y), -sin(rot.y), 0., sin(rot.y), cos(rot.y)) * pos; }

  void main() {
    float h_off = .6;
    vec3 p = rotate((pos-vec3(0, 0, h_off))*.9, mouse.xy*6.28);
    gl_Position = vec4(p, 1);
    vpos = pos; vcolor = color; vnormal = normal; vitr = itr;
  }`, `
  in float vitr;
  in vec3 vpos;
  in vec4 vcolor;
  in vec3 vnormal;
  out vec4 out_col;
  void main() {
    vec3 col = vcolor.xyz + vec3(0, 0, vcolor.w)*3.0;

    //vec3 light_dir = normalize(vec3(sin(vitr/100.), cos(vitr/100.), .5) - vpos);
    vec3 light_dir = normalize(vec3(-1, -1, 1) - vpos);
    float shading = max(0., dot(light_dir, vnormal)) + .09;
    out_col = vec4(col * shading, 1);
  }`);

const tris_buf = gl.createBuffer();
const colors_buf = gl.createBuffer();
const normals_buf = gl.createBuffer();

const dt = 1.2;
const init_volume = 1.;
const evap_rate = .001;
const depo_rate = .9;
const friction = .9999;
const stream_mix = .01;

const params = [ dt, init_volume, evap_rate, depo_rate, friction, stream_mix ];

for (let i = 0; i < 2000; i++) {
  wasm.update_physics(...params);
}

let itr = 0; const loop = () => { itr += 1; window.requestAnimationFrame(loop);

  for (let i = 0; i < 100; i++) wasm.update_physics(...params);
  if (itr % 10 == 0) {
    console.log(2000 + itr*100);

    wasm.update_geometry();
    gpu.upload_arr(gl, tris, tris_buf, 0);
    gpu.upload_arr(gl, colors, colors_buf, 1, 4);
    gpu.upload_arr(gl, normals, normals_buf, 2);
  }

  gl.useProgram(draw);
  gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
  gl.uniform3f(gl.getUniformLocation(draw, 'mouse'), mx, my, mouse_active);
  gl.uniform1f(gl.getUniformLocation(draw, 'itr'), itr);
  gl.drawArrays(gl.TRIANGLES, 0, tris.length/3);

[mvx, mvy] = [mvx*.98, mvy*.98]}; loop();

