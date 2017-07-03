uniform mat4 matrix;

in vec2 v;
in vec2 tc;
out vec2 v_texcoord;

void main() {
  gl_Position = matrix * vec4(v.x, v.y, 0.0, 1.0);
  v_texcoord = tc;
}
