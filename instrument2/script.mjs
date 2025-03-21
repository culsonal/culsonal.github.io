
import * as gpu from './gpu.mjs';

let actx = null;
window.document.onpointerdown = e => {
  if (actx == null) actx = new AudioContext({ latencyHint: 'interactive' }); }
const play_buf = (buf, len=-1) => {
  if (actx == null) return;
  const abuf = actx.createBuffer(1, len < 0 ? buf.length : len, 44100);
  abuf.copyToChannel(buf, 0);
  const src = actx.createBufferSource();
  src.buffer = abuf;
  src.connect(actx.destination);
  src.start();
};

const gl = document.getElementById('c').getContext('webgl2');
gl.getExtension("EXT_color_buffer_float");
gl.getExtension("EXT_color_buffer_half_float");
gl.canvas.style.width  = Math.min(window.innerWidth, window.innerHeight);
gl.canvas.style.height = Math.min(window.innerWidth, window.innerHeight);
gl.canvas.width = 210//parseInt(gl.canvas.style.width)/3;
gl.canvas.height = 210//parseInt(gl.canvas.style.height)/3;
console.log(gl.canvas.width, gl.canvas.height);

const update = gpu.make_frag_prog(gl, `
  out vec4 out_col;
  #define DIM (ivec2(${gl.canvas.width}, ${gl.canvas.height}))
  ${await (await fetch('sdf.glsl')).text()}
  ${await (await fetch('hash.glsl')).text()}
  ${await (await fetch('shader.glsl')).text()}
  void main() { update(); }`);

const draw = gpu.make_frag_prog(gl, `
  out vec4 out_col;
  #define DIM (ivec2(${gl.canvas.width}, ${gl.canvas.height}))
  ${await (await fetch('sdf.glsl')).text()}
  ${await (await fetch('hash.glsl')).text()}
  ${await (await fetch('shader.glsl')).text()}
  void main() { draw(); }`);

const record = gpu.make_frag_prog(gl, `
  #define RECORD
  out float out_col;
  #define DIM (ivec2(${gl.canvas.width}, ${gl.canvas.height}))
  ${await (await fetch('sdf.glsl')).text()}
  ${await (await fetch('hash.glsl')).text()}
  ${await (await fetch('shader.glsl')).text()}
  void main() { record(); }`);


const rec_dim = [ gl.canvas.width, gl.canvas.height ];
let [ tex1, fb1 ] = gpu.make_tex_and_fb(gl, [ gl.canvas.width, gl.canvas.height ]);
let [ tex2, fb2 ] = gpu.make_tex_and_fb(gl, [ gl.canvas.width, gl.canvas.height ]);
let [ rec_tex, rec_fb ] = gpu.make_tex_and_fb(gl, rec_dim, gl.R32F, gl.RED);
const rec_pbo = gl.createBuffer();

let [ mx, my, mpx, mpy, mvx, mvy, mouse_active ] = [ .5, .5, .5, .5, 0, 0, 0 ];
gl.canvas.onpointerdown = e => mouse_active = 1;
gl.canvas.onpointerup   = e => mouse_active = 0;
gl.canvas.onpointermove = e => { e.preventDefault();
  [ mpx, mpy, mx, my ] = [ mx, my,
    e.offsetX/parseInt(gl.canvas.style.width), 1-e.offsetY/parseInt(gl.canvas.style.height) ];
  [ mvx, mvy ] = [ mx-mpx, my-mpy ]; };

const run_sim = (num_sim_steps=1000, max_num_sim_steps=10_000, start_itr=0) => {
  let itr = start_itr;
  for (let i = 0; i < num_sim_steps; i++) {
    gpu.capture_with_framebuffer(gl, fb2, [ gl.canvas.width, gl.canvas.height ], () => {
      gl.useProgram(update);
      gpu.bind_tex(gl, update, tex1, 'img');
      gl.uniform1i(gl.getUniformLocation(update, 'num_sim_steps'), max_num_sim_steps);
      gl.uniform3f(gl.getUniformLocation(update, 'mouse'), mx, my, mouse_active);
      gl.uniform2f(gl.getUniformLocation(update, 'mouse_vel'), mvx, mvy);
      gl.uniform1i(gl.getUniformLocation(update, 'itr'), itr);
      gl.drawArrays(gl.TRIANGLES, 0, 6);
    });

    [ fb1, tex1, fb2, tex2 ] = [ fb2, tex2, fb1, tex1 ];

    itr += 1;
  }

  return itr;
};

const play_recording = (num_sim_steps) => {
  gpu.capture_with_framebuffer(gl, rec_fb, rec_dim, () => {
    gl.useProgram(record);
    gpu.bind_tex(gl, record, tex1, 'img');
    gl.uniform1i(gl.getUniformLocation(record, 'num_sim_steps'), num_sim_steps);
    gl.uniform3f(gl.getUniformLocation(record, 'mouse'), mx, my, mouse_active);
    gl.uniform2f(gl.getUniformLocation(record, 'mouse_vel'), mvx, mvy);
    gl.uniform1i(gl.getUniformLocation(record, 'itr'), itr);
    gl.drawArrays(gl.TRIANGLES, 0, 6);
  });

  gl.bindFramebuffer(gl.FRAMEBUFFER, rec_fb);
  const rec_rows = Math.ceil(num_sim_steps / gl.canvas.width);
  const data = new Float32Array(gl.canvas.width*rec_rows);
  gl.readPixels(
    0, 0,
    gl.canvas.width,
    rec_rows,
    gl.RED, gl.FLOAT,
    data
  );
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);

  play_buf(data, num_sim_steps);
  //console.log(data);
};

// could be used for auto-tuning later
const run_benchmark = () => {

  let log = '';
  const num_steps = [ 5, 6, 7, 8, 9, 10, 50, 500, 1000, 2000, 5000, 10_000, 20_000 ];
  let best_opt = [ 10e10, -1 ];
  for (let i = 0; i < num_steps.length; i++) {
    gl.flush();

    const t1 = performance.now();
    run_sim(num_steps[i], 0);
    const t2 = performance.now();
    const works = (t2-t1)/1000 < num_steps[i]/44100 ? '✅' : '❌';
    log += `${works} (${gl.canvas.width}x${gl.canvas.height}) chunking size ${num_steps[i]}:`+
      `took ${(t2-t1)/1000} seconds to run ${num_steps[i]/44100} seconds of simulation \n`;

    const score = (t2-t1)/1000 - num_steps[i]/44100;
    if (score < best_opt[0]) best_opt = [ score, num_steps[i] ];
  }

  console.log(log);
  console.log(`best option: ${best_opt[1]} with ${best_opt[0]} difference.`);
};

//run_benchmark();

let itr = 0; const loop = () => { window.requestAnimationFrame(loop);

  // physics
  const num_sim_steps = 500;
  const max_num_sim_steps = 44100/2;
  itr = run_sim(num_sim_steps, max_num_sim_steps, itr);

  gl.useProgram(draw);
  gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
  gpu.bind_tex(gl, draw, tex1, 'img');
  gl.uniform1i(gl.getUniformLocation(draw, 'num_sim_steps'), max_num_sim_steps);
  gl.uniform3f(gl.getUniformLocation(draw, 'mouse'), mx, my, mouse_active);
  gl.uniform2f(gl.getUniformLocation(draw, 'mouse_vel'), mvx, mvy);
  gl.uniform1f(gl.getUniformLocation(draw, 'itr'), itr);
  gl.drawArrays(gl.TRIANGLES, 0, 6);

  if (itr > max_num_sim_steps) {
    console.log('playing audio');
    play_recording(max_num_sim_steps);
    itr = 1;
  }

[mvx, mvy] = [mvx*.98, mvy*.98]}; loop();


