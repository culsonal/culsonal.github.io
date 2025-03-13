
import * as gpu from './gpu.mjs';

const gl = document.getElementById('c').getContext('webgl2');
gl.getExtension("EXT_color_buffer_float");
gl.canvas.style.width = window.innerWidth;
gl.canvas.style.height = window.innerHeight;
gl.canvas.width = parseInt(gl.canvas.style.width)/4;
gl.canvas.height = parseInt(gl.canvas.style.height)/4;

let [ tex1, fb1 ] = gpu.make_tex_and_fb(gl, [ gl.canvas.width, gl.canvas.height ]);
let [ tex2, fb2 ] = gpu.make_tex_and_fb(gl, [ gl.canvas.width, gl.canvas.height ]);

const update = gpu.make_frag_prog(gl, `
  #define DIM (ivec2(${gl.canvas.width}, ${gl.canvas.height}))
  ${await (await fetch('hash.glsl')).text()}
  ${await (await fetch('reintegration_tracking.glsl')).text()}
  void main() { update(); }`);

const draw = gpu.make_frag_prog(gl, `
  #define DIM (ivec2(${gl.canvas.width}, ${gl.canvas.height}))
  ${await (await fetch('hash.glsl')).text()}
  ${await (await fetch('reintegration_tracking.glsl')).text()}
  void main() { draw(); }`);

//let [ mx, my, mpx, mpy, mvx, mvy ] = [ .5, .5, .5, .5, 0, 0 ];
let [ mx, my, mpx, mpy, mvx, mvy ] = [ 0, 0, 0, 0, 0, 0 ];
gl.canvas.onpointermove = e => { e.preventDefault();
  [ mpx, mpy, mx, my ] = [ mx, my,
    e.offsetX/parseInt(gl.canvas.style.width), 1-e.offsetY/parseInt(gl.canvas.style.height) ];
  [ mvx, mvy ] = [ mx-mpx, my-mpy ]; };

let itr = -1; const loop = () => { itr += 1; window.requestAnimationFrame(loop);

for (let i = 0; i < 12; i++) {
  // update physics
  gpu.capture_with_framebuffer(gl, fb2, [ gl.canvas.width, gl.canvas.height ], () => {
    gl.useProgram(update);
    gpu.bind_tex(gl, update, tex1, 'img');
    gl.uniform4f(gl.getUniformLocation(update, 'mouse'), mx, my, mvx, mvy);
    gl.uniform1i(gl.getUniformLocation(update, 'itr'), itr);
    gl.drawArrays(gl.TRIANGLES, 0, 6);
  });

  [ fb1, tex1, fb2, tex2 ] = [ fb2, tex2, fb1, tex1 ];

  itr += 1;
}

// just for seeing state
gl.useProgram(draw);
gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
gpu.bind_tex(gl, draw, tex1, 'img');
gl.uniform4f(gl.getUniformLocation(draw, 'mouse'), mx, my, mvx, mvy);
gl.uniform1i(gl.getUniformLocation(draw, 'itr'), itr);
gl.drawArrays(gl.TRIANGLES, 0, 6);


[mvx, mvy] = [mvx*.98, mvy*.98]}; loop();


