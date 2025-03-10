
export const make_prog_basic = (gl, vsrc, fsrc) => {
  const [ p, vs, fs ] = [ gl.createProgram(), gl.createShader(gl.VERTEX_SHADER), gl.createShader(gl.FRAGMENT_SHADER) ];
  gl.shaderSource(vs, vsrc); gl.compileShader(vs); gl.attachShader(p, vs);
  gl.shaderSource(fs, fsrc); gl.compileShader(fs); gl.attachShader(p, fs);
  gl.linkProgram(p); return p;
};

export const make_prog = (gl, vsrc, fsrc, precision='highp') => {
  return make_prog_basic(gl,
    `#version 300 es \n precision ${precision} float; precision ${precision} int;`+vsrc,
    `#version 300 es \n precision ${precision} float; precision ${precision} int;`+fsrc);
};

export const make_frag_prog = (gl, fsrc, precision='highp') => {
  return make_prog(gl, `void main() {
    if (gl_VertexID == 0) gl_Position = vec4(-1, -1, 0, 1);
    if (gl_VertexID == 1) gl_Position = vec4( 1, -1, 0, 1);
    if (gl_VertexID == 2) gl_Position = vec4(-1,  1, 0, 1);
    if (gl_VertexID == 3) gl_Position = vec4( 1,  1, 0, 1);
    if (gl_VertexID == 4) gl_Position = vec4(-1,  1, 0, 1);
    if (gl_VertexID == 5) gl_Position = vec4( 1, -1, 0, 1);
  }`, fsrc, precision);
};

export const make_tex = (gl, dims, ifrmt=gl.RGBA32F, frmt=gl.RGBA, dtype=gl.FLOAT) => {
  const type = dims.length == 2 ? gl.TEXTURE_2D : gl.TEXTURE_3D;
  const tex_buf = gl.createTexture();
  gl.bindTexture(type, tex_buf);
  if (dims.length == 2) gl.texImage2D(type, 0, ifrmt, dims[0], dims[1],          0, frmt, dtype, null);
  if (dims.length == 3) gl.texImage3D(type, 0, ifrmt, dims[0], dims[1], dims[2], 0, frmt, dtype, null);
  gl.texParameteri(type, gl.TEXTURE_MIN_FILTER, gl.NEAREST); gl.texParameteri(type, gl.TEXTURE_WRAP_S, gl.REPEAT);
  gl.texParameteri(type, gl.TEXTURE_MAG_FILTER, gl.NEAREST); gl.texParameteri(type, gl.TEXTURE_WRAP_T, gl.REPEAT);
  gl.bindTexture(type, null);
  return tex_buf;
};

export const make_fb = (gl, tex) => {
  const fb = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
  gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  return fb;
};

export const make_tex_and_fb = (gl, dims, ifrmt=gl.RGBA32F, frmt=gl.RGBA, dtype=gl.FLOAT) => {
  const tex = make_tex(gl, dims, ifrmt, frmt, dtype);
  return [ tex, make_fb(gl, tex) ];
};

export const bind_tex = (gl, prog, tex, name, idx=0) => {
  gl.uniform1i(gl.getUniformLocation(prog, name), idx);
  gl.activeTexture(gl.TEXTURE0 + idx);
  gl.bindTexture(gl.TEXTURE_2D, tex);
};

export const capture_with_framebuffer = (gl, fb, [ width, height ], callback) => {
  gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
  gl.viewport(0, 0, width, height);
  callback();
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
};

