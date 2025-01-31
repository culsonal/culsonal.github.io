
const draw_bodies = (ctx, [ bodies, num_bodies, num_pts_per_body, num_species ]) => {
  ctx.fillStyle = '#fff'; ctx.strokeStyle = '#fff';
  for (let i = 0; i < num_bodies; i++) {
    //ctx.strokeStyle = i % 2 == 0 ? '#0f0' : '#00f';
    //ctx.fillStyle = i % 2 == 0 ? '#0f0' : '#00f';

    ctx.beginPath();
    for (let j = 0; j < num_pts_per_body; j++) {
      const ii = i*num_pts_per_body + j;
      const [ x, y ] = bodies.subarray(ii*2*4, ii*2*4+2);
      ctx.lineTo((x*.5+.5)*ctx.canvas.width, (y*.5+.5)*ctx.canvas.height);
      if (j == 0) ctx.fillRect((x*.5+.5)*ctx.canvas.width-1.5, (y*.5+.5)*ctx.canvas.height-1.5, 3, 3);
    }
    ctx.stroke();
  }
};

const draw_map_v1 = (ctx, [ map, map_dim, num_species ], viz=-1) => {
  const num_map_vals = num_species+1;
  const img_data = ctx.getImageData(0, 0, map_dim, map_dim);
  for (let i = 0; i < map_dim; i++) for (let j = 0; j < map_dim; j++) {
    const map_vals = map.subarray((i*map_dim + j)*num_map_vals, (i*map_dim + j)*num_map_vals+num_map_vals);
    const [ r, g, b ] = map_vals.subarray(0, 3); // assuming 2 species and 1 nutrient thing for now
    let pix;
    if (viz == -1) pix = [ r, g, b, 255 ];
    else {
      const v = map_vals[viz]*255;
      //pix = [ v, v, v, 255 ];
      pix = [ 0, v, 0, 255 ];
    }
    img_data.data.set(pix, (i*map_dim + j)*4);
  }
  ctx.putImageData(img_data, 0, 0);
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(ctx.canvas, 0, 0, map_dim, map_dim, 0, 0, ctx.canvas.width, ctx.canvas.height);
};

const draw_map_v2 = (ctx, [ map, map_dim, num_species ]) => {
  const num_map_vals = num_species+1;
  const img_data1 = ctx.getImageData(0, 0, map_dim, map_dim);
  const img_data2 = ctx.getImageData(map_dim,   0, map_dim, map_dim);
  const img_data3 = ctx.getImageData(map_dim*2, 0, map_dim, map_dim);
  for (let i = 0; i < map_dim; i++) for (let j = 0; j < map_dim; j++) {
    const map_vals = map.subarray((i*map_dim + j)*num_map_vals, (i*map_dim + j)*num_map_vals+num_map_vals);
    const [ r, g, b ] = map_vals.subarray(0, 3); // assuming 2 species and 1 nutrient thing for now
    img_data1.data.set([ r*255, 0, 0, 255 ], (i*map_dim + j)*4);
    img_data2.data.set([ 0, g*255, 0, 255 ], (i*map_dim + j)*4);
    img_data3.data.set([ 0, 0, b*255, 255 ], (i*map_dim + j)*4);
  }
  ctx.putImageData(img_data1, 0, 0);
  ctx.putImageData(img_data2, map_dim, 0);
  ctx.putImageData(img_data3, map_dim*2, 0);
};

const wasm_ = await WebAssembly.instantiateStreaming(fetch('a.out.wasm'));
const wasm = wasm_.instance.exports;
wasm.init(Math.floor(Math.random()*100), .2);

const num_bodies = wasm.get_num_bodies();
const num_pts_per_body = wasm.get_num_pts_per_body();
const bodies_len = num_bodies * num_pts_per_body * 2 * 4;
const bodies = new Float32Array(wasm.memory.buffer, wasm.get_bodies_ptr(), bodies_len);
const map_dim = wasm.get_map_dim();
const num_species = wasm.get_num_species();
const map = new Float32Array(wasm.memory.buffer, wasm.get_map_ptr(), map_dim*map_dim*(num_species+1));

const ctx = document.getElementById('c').getContext('2d', { willReadFrequently: true });
ctx.canvas.width = ctx.canvas.height = Math.min(window.innerWidth, window.innerHeight)-2;

const viz = location.hash.length == 0 ? 0 : parseInt(location.hash.slice(1));
const loop = () => { window.requestAnimationFrame(loop);

  ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
  //draw_map_v1(ctx, [ map, map_dim, num_species ], 1);
  //draw_map_v1(ctx, [ map, map_dim, num_species ], viz);
  //draw_map_v2(ctx, [ map, map_dim, num_species ]);
  draw_bodies(ctx, [ bodies, num_bodies, num_pts_per_body, num_species ]);
  wasm.update();

}; loop();


