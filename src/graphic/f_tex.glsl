uniform vec4 colour;
uniform sampler2D tex;

varying vec2 v_texcoord;

void main() {
  gl_FragColor = texture2D(tex, v_texcoord);
}
