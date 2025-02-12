
// Maybe merge this into the main lib eventually

export const good_defaults_gl = (canvas, depth=true, full_screen=true, canvas_div=1) => {
  if (full_screen) {
    canvas.width  = Math.min(window.innerWidth, window.innerHeight)/canvas_div;
    canvas.height = Math.min(window.innerWidth, window.innerHeight)/canvas_div;
    canvas.style = window.innerHeight < window.innerWidth ? 'height:100%;width:auto' : 'width:100%;height:auto';
  }
  const gl = canvas.getContext('webgl2', {
    alpha: false, desynchronized: true, antialias: true,
    powerPreference: 'high-performance', //preserveDrawingBuffer: true,
  }); 
  gl.viewport(0, 0, canvas.width, canvas.height);
  if (depth) gl.enable(gl.DEPTH_TEST);
  gl.getExtension("EXT_color_buffer_float"); gl.getExtension("EXT_color_buffer_half_float");
  //gl.enable(gl.CULL_FACE); gl.cullFace(gl.FRONT_AND_BACK); gl.cullFace(gl.FRONT); gl.cullFace(gl.BACK); 
  return gl;
};

export const make_prog_basic = (gl, vsrc, fsrc, tf_out=[], tf_type=gl.SEPARATE_ATTRIBS) => {
  const [ p, vs, fs ] = [ gl.createProgram(), gl.createShader(gl.VERTEX_SHADER), gl.createShader(gl.FRAGMENT_SHADER) ];
  gl.shaderSource(vs, vsrc); gl.compileShader(vs); gl.attachShader(p, vs);
  gl.shaderSource(fs, fsrc); gl.compileShader(fs); gl.attachShader(p, fs);
  if (tf_out.length > 0) gl.transformFeedbackVaryings(p, tf_out, tf_type);
  gl.linkProgram(p); return p;
};

export const make_prog = (gl, vsrc, fsrc='', precision='highp', tf_out=[], tf_type=gl.SEPARATE_ATTRIBS) => {
  return make_prog_basic(gl,
    `#version 300 es \n precision ${precision} float; precision ${precision} int;`+vsrc,
    `#version 300 es \n precision ${precision} float; precision ${precision} int;`+(fsrc=='' ? 'void main(){}' : fsrc),
    tf_out, tf_type);
};

export const make_frag_prog = (gl, tex_dim, fsrc, precision='highp') => {
  return make_prog(gl,
    `void main() { gl_Position = vec4(0, 0, 0, 1); gl_PointSize = float(${tex_dim}); }`,
    fsrc, precision);
};

export const make_tex = (gl, dims, ifrmt=gl.RGBA32F, frmt=gl.RGBA, dtype=gl.FLOAT) => {
  const type = dims.length == 2 ? gl.TEXTURE_2D : gl.TEXTURE_3D;
  const tex_buf = gl.createTexture();
  gl.bindTexture(type, tex_buf);
  if (dims.length == 2) gl.texImage2D(type, 0, ifrmt, dims[0], dims[1],          0, frmt, dtype, null);
  if (dims.length == 3) gl.texImage3D(type, 0, ifrmt, dims[0], dims[1], dims[2], 0, frmt, dtype, null);
  gl.texParameteri(type, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(type, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texParameteri(type, gl.TEXTURE_WRAP_S, gl.REPEAT);
  gl.texParameteri(type, gl.TEXTURE_WRAP_T, gl.REPEAT);
  //gl.texParameteri(type, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  //gl.texParameteri(type, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  gl.bindTexture(type, null);
  return tex_buf;
};

export const make_framebuffer = (gl, tex_bufs) => {
  const fb = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
  for (let i = 0; i < tex_bufs.length; i++)
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0+i, gl.TEXTURE_2D, tex_bufs[i], 0);
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  return fb;
};

export const make_tex2 = (gl, dims, ifrmt=gl.RGBA32F, frmt=gl.RGBA, dtype=gl.FLOAT) => {
  return {
    A: make_tex(gl, dims, ifrmt, frmt, dtype),
    B: make_tex(gl, dims, ifrmt, frmt, dtype)
  };
};

export const make_framebuffer2 = (gl, tex_bufs) => {
  return {
    A: make_framebuffer(gl, tex_bufs.map(el => el.A)),
    B: make_framebuffer(gl, tex_bufs.map(el => el.B))
  };
};

export const swap_framebuffer = (gl, fb, tex) => {
  [ fb.A, fb.B, tex.A, tex.B ] = [ fb.B, fb.A, tex.B, tex.A ];
};

export const make_array_buf = (gl, arr, hint=gl.DYNAMIC_DRAW) => {
  const buf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  gl.bufferData(gl.ARRAY_BUFFER, arr, hint);
  gl.bindBuffer(gl.ARRAY_BUFFER, null);
  return buf;
};

export const make_empty_array_buf = (gl, size, bytes_per_elem=4, hint=gl.DYNAMIC_DRAW) => {
  const buf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  gl.bufferData(gl.ARRAY_BUFFER, size*bytes_per_elem, hint);
  gl.bindBuffer(gl.ARRAY_BUFFER, null);
  return buf;
};

export const make_pbo = (gl, size, bytes_per_elem=4, hint=gl.DYNAMIC_DRAW) => {
  const pbo = gl.createBuffer();
  gl.bindBuffer(gl.PIXEL_UNPACK_BUFFER, pbo);
  gl.bufferData(gl.PIXEL_UNPACK_BUFFER, size*bytes_per_elem, gl.DYNAMIC_DRAW);
  gl.bindBuffer(gl.PIXEL_UNPACK_BUFFER, null);
  return pbo;
};

export const copy_buf_to_tex_via_pbo = (gl, buf, pbo, tex, [ width, height, internal_format, format, dtype ],
    size, bytes_per_float=4) => {
  gl.bindBuffer(gl.COPY_READ_BUFFER, buf);
  gl.bindBuffer(gl.PIXEL_UNPACK_BUFFER, pbo);
  gl.copyBufferSubData(gl.COPY_READ_BUFFER, gl.PIXEL_UNPACK_BUFFER, 0, 0, bytes_per_float*size);
  gl.bindBuffer(gl.PIXEL_UNPACK_BUFFER, null);
  gl.bindBuffer(gl.COPY_READ_BUFFER, null);

  gl.bindTexture(gl.TEXTURE_2D, tex);
  gl.bindBuffer(gl.PIXEL_UNPACK_BUFFER, pbo);
  gl.texImage2D(gl.TEXTURE_2D, 0, internal_format, width, height, 0, format, dtype, 0);
  gl.bindBuffer(gl.PIXEL_UNPACK_BUFFER, null);
  gl.bindTexture(gl.TEXTURE_2D, null);
};

// stuff you'd use for multiple outputs
//  gl.framebufferTextureLayer(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0+i, tex, 0, i);
//  gl.drawBuffers([ gl.COLOR_ATTACHMENT0, gl.COLOR_ATTACHMENT1, ... ]);
export const capture_with_framebuffer = (gl, fb, [ width, height ], callback) => {
  gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
  gl.viewport(0, 0, width, height);
  //gl.clear(gl.COLOR_BUFFER_BIT); // or make not default?
  callback();
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
};

// missing: element array buffer? or do that inside 'draw_indexed_arrays' maybe
export const setup_prog_inputs_and_outputs = (gl, prog, inputs={}, outputs={}) => {
  const num_inputs = gl.getProgramParameter(prog, gl.ACTIVE_ATTRIBUTES);
  for (let i = 0; i < num_inputs; i++) { // missing: bool type, matrix types
    const input = gl.getActiveAttrib(prog, i);
    if (input.name == 'gl_VertexID') continue; // apparantely this comes in as an attribute, weird
    const t = input.type; let [ dtype, dim ] = [ 0, 0 ];
    if (t == gl.UNSIGNED_INT)      [ dtype, dim ] = [ gl.UNSIGNED_INT, 1 ];
    if (t == gl.UNSIGNED_INT_VEC2) [ dtype, dim ] = [ gl.UNSIGNED_INT, 2 ];
    if (t == gl.UNSIGNED_INT_VEC3) [ dtype, dim ] = [ gl.UNSIGNED_INT, 3 ];
    if (t == gl.UNSIGNED_INT_VEC4) [ dtype, dim ] = [ gl.UNSIGNED_INT, 4 ];
    if (t == gl.FLOAT)      [ dtype, dim ] = [ gl.FLOAT, 1 ];
    if (t == gl.FLOAT_VEC2) [ dtype, dim ] = [ gl.FLOAT, 2 ];
    if (t == gl.FLOAT_VEC3) [ dtype, dim ] = [ gl.FLOAT, 3 ];
    if (t == gl.FLOAT_VEC4) [ dtype, dim ] = [ gl.FLOAT, 4 ];
    if (t == gl.INT)        [ dtype, dim ] = [ gl.INT, 1 ];
    if (t == gl.INT_VEC2)   [ dtype, dim ] = [ gl.INT, 2 ];
    if (t == gl.INT_VEC3)   [ dtype, dim ] = [ gl.INT, 3 ];
    if (t == gl.INT_VEC4)   [ dtype, dim ] = [ gl.INT, 4 ];
    const loc = gl.getAttribLocation(prog, input.name);
    gl.bindBuffer(gl.ARRAY_BUFFER, inputs[input.name]);
    gl.enableVertexAttribArray(loc);
    gl.vertexAttribPointer(loc, dim, dtype, false, 0, 0);
  }

  let sampler_idx = 0;
  const num_uniforms = gl.getProgramParameter(prog, gl.ACTIVE_UNIFORMS);
  for (let i = 0; i < num_uniforms; i++) {
    const uniform = gl.getActiveUniform(prog, i);
    const loc = gl.getUniformLocation(prog, uniform.name);
    const [ t, d ] = [ uniform.type, inputs[uniform.name] ];
    // missing: FLOAT_MAT2, FLOAT_MAT3, FLOAT_MAT4, SAMPLER_CUBE
    if (t == gl.UNSIGNED_INT)      gl.uniform1ui(loc, d);
    if (t == gl.UNSIGNED_INT_VEC2) gl.uniform2ui(loc, d[0], d[1]);
    if (t == gl.UNSIGNED_INT_VEC3) gl.uniform3ui(loc, d[0], d[1], d[2]);
    if (t == gl.UNSIGNED_INT_VEC4) gl.uniform4ui(loc, d[0], d[1], d[2], d[3]);
    if (t == gl.FLOAT)      gl.uniform1f(loc, d);
    if (t == gl.FLOAT_VEC2) gl.uniform2f(loc, d[0], d[1]);
    if (t == gl.FLOAT_VEC3) gl.uniform3f(loc, d[0], d[1], d[2]);
    if (t == gl.FLOAT_VEC4) gl.uniform4f(loc, d[0], d[1], d[2], d[3]);
    if (t == gl.INT)        gl.uniform1i(loc, d);
    if (t == gl.INT_VEC2)   gl.uniform2i(loc, d[0], d[1]);
    if (t == gl.INT_VEC3)   gl.uniform3i(loc, d[0], d[1], d[2]);
    if (t == gl.INT_VEC4)   gl.uniform4i(loc, d[0], d[1], d[2], d[3]);
    if (t == gl.BOOL)       gl.uniform1b(loc, d); // just guessing it's uniformNb
    if (t == gl.BOOL_VEC2)  gl.uniform2b(loc, d[0], d[1]);
    if (t == gl.BOOL_VEC3)  gl.uniform3b(loc, d[0], d[1], d[2]);
    if (t == gl.BOOL_VEC4)  gl.uniform4b(loc, d[0], d[1], d[2], d[3]);
    if (t == gl.SAMPLER_2D || t == gl.INT_SAMPLER_2D || t == gl.UNSIGNED_INT_SAMPLER_2D ||
        t == gl.SAMPLER_3D || t == gl.INT_SAMPLER_3D || t == gl.UNSIGNED_INT_SAMPLER_3D) {
      let tex_type;
      if (t == gl.SAMPLER_2D || t == gl.INT_SAMPLER_2D || t == gl.UNSIGNED_INT_SAMPLER_2D) tex_type = gl.TEXTURE_2D;
      if (t == gl.SAMPLER_3D || t == gl.INT_SAMPLER_3D || t == gl.UNSIGNED_INT_SAMPLER_3D) tex_type = gl.TEXTURE_3D;
      gl.uniform1i(loc, sampler_idx);
      gl.activeTexture(gl.TEXTURE0 + sampler_idx);
      gl.bindTexture(tex_type, d);
      sampler_idx += 1;
    }
  }

  const num_outputs = gl.getProgramParameter(prog, gl.TRANSFORM_FEEDBACK_VARYINGS);
  for (let i = 0; i < num_outputs; i++) {
    const tf_name = gl.getTransformFeedbackVarying(prog, i).name;
    gl.bindBufferBase(gl.TRANSFORM_FEEDBACK_BUFFER, i, outputs[tf_name]);
  }
};

export const run_prog_pts = (gl, prog, num_pts, inputs={}, outputs={}, dont_render_if_tf=true) => {
  const is_tf_prog = Object.keys(outputs).length > 0;
  gl.useProgram(prog);
  setup_prog_inputs_and_outputs(gl, prog, inputs, outputs);
  if (is_tf_prog) { gl.beginTransformFeedback(gl.POINTS); if (dont_render_if_tf) gl.enable(gl.RASTERIZER_DISCARD); }

  gl.drawArrays(gl.POINTS, 0, num_pts);

  if (is_tf_prog) { gl.endTransformFeedback(); if (dont_render_if_tf) gl.disable(gl.RASTERIZER_DISCARD); }
  for (let i = 0; i < Object.keys(outputs).length; i++) gl.bindBufferBase(gl.TRANSFORM_FEEDBACK_BUFFER, i, null);
  gl.bindBuffer(gl.ARRAY_BUFFER, null);
};

