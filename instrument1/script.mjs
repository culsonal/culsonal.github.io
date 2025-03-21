
import * as gpu from './gpu.mjs';

const play_buf = buf => {
  const actx = new AudioContext();
  const abuf = new AudioBuffer({ length: buf.length, sampleRate: 44100 });
  abuf.getChannelData(0).set(buf);
  const src = actx.createBufferSource();
  src.buffer = abuf;
  src.connect(actx.destination);
  src.start();
};

const gl = document.getElementById('c').getContext('webgl2');
gl.getExtension("EXT_color_buffer_float");
gl.canvas.style.width  = Math.min(window.innerWidth, window.innerHeight);
gl.canvas.style.height = Math.min(window.innerWidth, window.innerHeight);
gl.canvas.width = 410//220//parseInt(gl.canvas.style.width)/3;
gl.canvas.height = 410//110//parseInt(gl.canvas.style.height)/3;

const update = gpu.make_frag_prog(gl, `
  out vec4 out_col;
  #define DIM (ivec2(${gl.canvas.width}, ${gl.canvas.height}))
  ${await (await fetch('hash.glsl')).text()}
  ${await (await fetch('shader.glsl')).text()}
  void main() { update(); }`);

const draw = gpu.make_frag_prog(gl, `
  out vec4 out_col;
  #define DIM (ivec2(${gl.canvas.width}, ${gl.canvas.height}))
  ${await (await fetch('hash.glsl')).text()}
  ${await (await fetch('shader.glsl')).text()}
  void main() { draw(); }`);

const record = gpu.make_frag_prog(gl, `
  #define RECORD
  out float out_col;
  #define DIM (ivec2(${gl.canvas.width}, ${gl.canvas.height}))
  ${await (await fetch('hash.glsl')).text()}
  ${await (await fetch('shader.glsl')).text()}
  void main() { record(); }`);

const gldim = [ gl.canvas.width, gl.canvas.height ];
let [ tex1, fb1 ] = gpu.make_tex_and_fb(gl, gldim);
let [ tex2, fb2 ] = gpu.make_tex_and_fb(gl, gldim);
let [ rec_tex, rec_fb ] = gpu.make_tex_and_fb(gl, gldim, gl.R32F, gl.RED);

const handle_mouse = e => { e.preventDefault();
  [ mpx, mpy, mx, my ] = [ mx, my,
    e.offsetX/parseInt(gl.canvas.style.width), 1-e.offsetY/parseInt(gl.canvas.style.height) ];
  [ mvx, mvy ] = [ mx-mpx, my-mpy ]; 

  itr = 0; //// **** ////
};

let [ mx, my, mpx, mpy, mvx, mvy ] = [ .8, .5, .8, .5, 0, 0 ];
gl.canvas.onpointermove = e => handle_mouse(e);
window.document.onpointerdown = e => {
  gpu.capture_with_framebuffer(gl, rec_fb, [ gl.canvas.width, gl.canvas.height ], () => {
    gl.useProgram(record);
    gpu.bind_tex(gl, record, tex1, 'img');
    gl.uniform4f(gl.getUniformLocation(record, 'mouse'), mx, my, mvx, mvy);
    gl.uniform1i(gl.getUniformLocation(record, 'itr'), itr);
    gl.drawArrays(gl.TRIANGLES, 0, 6);
  });

  gl.bindFramebuffer(gl.FRAMEBUFFER, rec_fb);
  const data = new Float32Array(gl.canvas.width*gl.canvas.height);
  gl.readPixels(
    0, 0,
    gl.canvas.width,
    gl.canvas.height,
    gl.RED, gl.FLOAT,
    data
  );

  play_buf(data);
}

let itr = 0; const loop = () => { window.requestAnimationFrame(loop);

  // physics
  for (let i = 0; i < 600; i++) {
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

  //if (itr == 500) console.log(tex1)

  // draw
  gl.useProgram(draw);
  gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
  gpu.bind_tex(gl, draw, tex1, 'img');
  gl.uniform4f(gl.getUniformLocation(draw, 'mouse'), mx, my, mvx, mvy);
  gl.uniform1i(gl.getUniformLocation(draw, 'itr'), itr);
  gl.drawArrays(gl.TRIANGLES, 0, 6);

  //if (itr > gl.canvas.width*gl.canvas.height) itr = 1;

[mvx, mvy] = [mvx*.98, mvy*.98]}; loop();


