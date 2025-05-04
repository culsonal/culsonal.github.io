
import * as gpu from './gpu.mjs';

const gl = document.getElementById('c').getContext('webgl2');
gl.getExtension("EXT_color_buffer_float");
gl.canvas.style.width  = Math.min(window.innerWidth, window.innerHeight)-2;//innerWidth
gl.canvas.style.height = Math.min(window.innerWidth, window.innerHeight)-2;//innerHeight
gl.canvas.width = 512;//parseInt(gl.canvas.style.width)/4;
gl.canvas.height = 512;//parseInt(gl.canvas.style.height)/4;
let [ mx, my, mvx, mvy, mouse_active ] = [ .5, .5, 0, 0, 0 ];
gl.canvas.onpointerdown = e => mouse_active = 1;
gl.canvas.onpointerup   = e => mouse_active = 0;
gl.canvas.onpointermove = e => [ mvx, mvy, mx, my ] = [
  e.movementX/parseInt(gl.canvas.style.width), e.movementY/parseInt(gl.canvas.style.height),
  e.offsetX/parseInt(gl.canvas.style.width), 1-e.offsetY/parseInt(gl.canvas.style.height) ];

const draw = gpu.make_frag_prog(gl, `
  #define DIM (ivec2(${gl.canvas.width}, ${gl.canvas.height}))
  ${await (await fetch('hash.glsl')).text()}
  ${await (await fetch('sdf2.glsl')).text()}
  ${await (await fetch('sdf3.glsl')).text()}
  ${await (await fetch('shader.glsl')).text()}
  void main() { draw(); }`);

const focal_length = 2; // todo make this an input (slider or something)

let itr = 0; const loop = () => { itr++; window.requestAnimationFrame(loop);

  gl.useProgram(draw);
  gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
  gl.uniform3f(gl.getUniformLocation(draw, 'mouse'), mx, my, mouse_active);
  gl.uniform1f(gl.getUniformLocation(draw, 'itr'), itr);
  gl.uniform1f(gl.getUniformLocation(draw, 'focal_length'), focal_length);
  gl.drawArrays(gl.TRIANGLES, 0, 6);

[mvx, mvy] = [mvx*.98, mvy*.98]}; loop();


