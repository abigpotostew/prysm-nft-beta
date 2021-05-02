#ifdef GL_ES
precision mediump float;
#endif


uniform float u_time;
uniform vec2 u_resolution;

float dot2( in vec2 v ) { return dot(v,v); }
float dot2( in vec3 v ) { return dot(v,v); }
float ndot( in vec2 a, in vec2 b ) { return a.x*b.x - a.y*b.y; }
float opSmoothUnion( float d1, float d2, float k )
{
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h);
}

mat4 Rot4X(float a ) {
    float c = cos( a );
    float s = sin( a );
    return mat4( 1, 0, 0, 0,
                 0, c,-s, 0,
                 0, s, c, 0,
                 0, 0, 0, 1 );
}


// Return 4x4 rotation Y matrix
// angle in radians
// ========================================
mat4 Rot4Y(float a ) {
    float c = cos( a );
    float s = sin( a );
    return mat4( c, 0, s, 0,
                 0, 1, 0, 0,
                -s, 0, c, 0,
                 0, 0, 0, 1 );
}

// Return 4x4 rotation Z matrix
// angle in radians
// ========================================
mat4 Rot4Z(float a ) {
    float c = cos( a );
    float s = sin( a );
    return mat4(
        c,-s, 0, 0,
        s, c, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
     );
}

// d = distance to move
// ========================================
vec3 opTrans( vec3 p, vec3 d ) {
    return p - d;
}

// if no support for GLSL 1.2+
//     #version 120
// ========================================
mat4 transposeM4(in mat4 m ) {
    vec4 r0 = m[0];
    vec4 r1 = m[1];
    vec4 r2 = m[2];
    vec4 r3 = m[3];

    mat4 t = mat4(
         vec4( r0.x, r1.x, r2.x, r3.x ),
         vec4( r0.y, r1.y, r2.y, r3.y ),
         vec4( r0.z, r1.z, r2.z, r3.z ),
         vec4( r0.w, r1.w, r2.w, r3.w )
    );
    return t;
}

// Translate is simply: p - d
// opTx will do transpose(m)
// p' = m*p
//    = [m0 m1 m2 m3 ][ p.x ]
//      [m4 m5 m6 m7 ][ p.y ]
//      [m8 m9 mA mB ][ p.z ]
//      [mC mD mE mF ][ 1.0 ]
// ========================================
mat4 Loc4( vec3 p ) {
    p *= -1.;
    return mat4(
        1,  0,  0,  p.x,
        0,  1,  0,  p.y,
        0,  0,  1,  p.z,
        0,  0,  0,  1
    );
}

// Note: m must already be inverted!
// TODO: invert(m) transpose(m)
// Op Rotation / Translation
// ========================================
vec3 opTx( vec3 p, mat4 m ) {   // BUG in iq's docs, should be q
    return (transposeM4(m)*vec4(p,1.0)).xyz;
}


float sdSphere( vec3 p, float s )
{
  return length(p)-s;
} 


float sdTriPrism( vec3 p, vec2 h )
{
    float k = sqrt(3.0);
    h.x *= 0.5*k;
    p.xy /= h.x;
    p.x = abs(p.x) - 1.0;
    p.y = p.y + 1.0/k;
    if( p.x+k*p.y>0.0 ) p.xy=vec2(p.x-k*p.y,-k*p.x-p.y)/2.0;
    p.x -= clamp( p.x, -2.0, 0.0 );
    float d1 = length(p.xy)*sign(-p.y)*h.x;
    float d2 = abs(p.z)-h.y;
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}


float sdOctogonPrism( in vec3 p, in float r, float h )
{
  vec3 k = vec3(-0.9238795325,   // sqrt(2+sqrt(2))/2 
                       0.3826834323,   // sqrt(2-sqrt(2))/2
                       0.4142135623 ); // sqrt(2)-1 
  // reflections
  p = abs(p);
  p.xy -= 2.0*min(dot(vec2( k.x,k.y),p.xy),0.0)*vec2( k.x,k.y);
  p.xy -= 2.0*min(dot(vec2(-k.x,k.y),p.xy),0.0)*vec2(-k.x,k.y);
  // polygon side
  p.xy -= vec2(clamp(p.x, -k.z*r, k.z*r), r);
  vec2 d = vec2( length(p.xy)*sign(p.y), p.z-h );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

// c is the sin/cos of the desired cone angle
float sdSolidAngle(vec3 pos, vec2 c, float ra)
{
    vec2 p = vec2( length(pos.xz), pos.y );
    float l = length(p) - ra;
    float m = length(p - c*clamp(dot(p,c),0.0,ra) );
    return max(l,m*sign(c.y*p.x-c.x*p.y));
}


// arbitrary orientation
float sdCylinder(vec3 p, vec3 a, vec3 b, float r)
{
    vec3 pa = p - a;
    vec3 ba = b - a;
    float baba = dot(ba,ba);
    float paba = dot(pa,ba);

    float x = length(pa*baba-ba*paba) - r*baba;
    float y = abs(paba-baba*0.5)-baba*0.5;
    float x2 = x*x;
    float y2 = y*y*baba;
    float d = (max(x,y)<0.0)?-min(x2,y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));
    return sign(d)*sqrt(abs(d))/baba;
}

float sdRoundCone(vec3 p, vec3 a, vec3 b, float r1, float r2)
{
    // sampling independent computations (only depend on shape)
    vec3  ba = b - a;
    float l2 = dot(ba,ba);
    float rr = r1 - r2;
    float a2 = l2 - rr*rr;
    float il2 = 1.0/l2;
    
    // sampling dependant computations
    vec3 pa = p - a;
    float y = dot(pa,ba);
    float z = y - l2;
    float x2 = dot2( pa*l2 - ba*y );
    float y2 = y*y*l2;
    float z2 = z*z*l2;

    // single square root!
    float k = sign(rr)*rr*rr*x2;
    if( sign(z)*a2*z2 > k ) return  sqrt(x2 + z2)        *il2 - r2;
    if( sign(y)*a2*y2 < k ) return  sqrt(x2 + y2)        *il2 - r1;
                            return (sqrt(x2*a2*il2)+y*rr)*il2 - r1;
}

float map(vec3 p)
{
    float d = 2.0;

    vec3 q = p; // q is tmp variable for changing position for the next op.
    vec3 spheresRandomVec = vec3(52.5126, 64.62744, 632.25);
    for (int i = 0; i < 5; i++) {
        float fi = float(i+10);
        float time = u_time * (fract(fi * 412.531 + 0.513) - 0.5) * 2.0;
        q = p + sin(time + fi * spheresRandomVec) * vec3(2.0, 2.0, 0.8);
        d = opSmoothUnion(
            sdSphere(q, mix(0.5, 1.0, fract(fi * 412.531 + 0.5124))),
            d,
            0.6
        );
    }

    vec3 otherRandomVec = vec3(52.5126, 64.62744, 632.25);
    int i=1;
    float fi = float(i);
    float time = u_time * (fract(fi * 412.531 + 0.513) - 0.5) * 2.0;
    d = opSmoothUnion(
        sdSphere(
            p + sin(time/3.0 + fi * otherRandomVec) * vec3(2.0, 2.0, 0.8), 
            mix(0.5, 1.0, fract(fi * 412.531 + 0.5124))),
        d,
        0.5
    );
    
    i=2;
    fi = float(i);
    d=opSmoothUnion(
            sdSphere(
                p+sin(time/2.0+ fi * otherRandomVec) * vec3(2.0, 2.0, 0.8),
                0.75),//mix(0.5, 1.0, fract(fi * 412.531 + 0.5124))),
                //(sin(time/2.0+fi *195.26402)+1.0)/4.0+0.5,
                //(sin(time/2.0+fi*10.0)+1.0)/4.0+0.5),
            d, 
            0.5);
    
    i=11;
    fi = float(i);
    d = opSmoothUnion(
        sdSphere(
            p + sin(time/3.0 + fi * otherRandomVec) * vec3(2.0, 2.0, 0.8), 
            mix(0.2, .7, fract(fi * 412.531 + 0.5124))),
        d,
        0.5
    );
     
    /*d = opSmoothUnion(
            sdOctogonPrism(
                p+sin(time/3.0+ fi * vec3(52.5126, 64.62744, 632.25)) * vec3(2.0, 2.0, 0.8),
                (sin(time/2.0+fi *195.26402)+1.0)/4.0+0.5,
                (sin(time/2.0+fi*10.0)+1.0)/4.0+0.5),
            d, 
            0.6);*/
    
    i=4;
    fi=float(i);
    // this one flies way out
    d = opSmoothUnion(
            sdRoundCone(
                p+sin(time/3.0+ fi * otherRandomVec) * vec3(5.0, 3.0, 0.8),
                sin(time/3.0+ fi * otherRandomVec) * vec3(2.0, 2.0, 0.8),
                vec3(1.0),//sin(time/3.0+ fi * vec3(52.5126, 64.62744, 632.25)) * vec3(2.0, 2.0, 0.8),
                .2,
                .5
                //(sin(time/2.0+fi *195.26402)+1.0)/4.0+0.5
                //(sin(time/2.0+fi*10.0)+1.0)/4.0+0.5
                ),
            d, 
            0.5);
    i=5;
    fi=float(i);
    d = opSmoothUnion(
            sdRoundCone(
                p+sin(time/3.0+ fi * otherRandomVec) * vec3(2.4, 3.99, 1.8),
                vec3(0.1,0.1,0.1),//sin(time/3.0+ fi * vec2(52.5126, 64.62744)) * vec2(1.0, 1.0),
                sin(time/3.0+ fi * otherRandomVec) * vec3(2.0, 2.0, 0.8),
                0.5,
                0.32
                //(sin(time/2.0+fi *195.26402)+1.0)/4.0+0.5
                //(sin(time/2.0+fi*10.0)+1.0)/4.0+0.5
                ),
            d, 
            0.5);


    float prysmSize = 1.0;
    q = p;
    q = opTx( q, Loc4( vec3(0.0,0.0,-prysmSize*1.0) ));
    q = opTx( q, Rot4X( u_time/3.0 ) ); 
    q = opTx(q, Rot4Y(u_time/5.0));
    // move the prysm closer
    float prysmGoopyness=.9;
    d = opSmoothUnion(sdTriPrism(q+vec3(0.0,0.0,1.0), vec2(prysmSize)), d, prysmGoopyness);
    
    
    return d;
} 

vec3 calcNormal(  vec3 p )
{
    float h = 1e-5; // or some other value
    vec2 k = vec2(1,-1);
    return normalize( k.xyy*map( p + k.xyy*h ) + 
                      k.yyx*map( p + k.yyx*h ) + 
                      k.yxy*map( p + k.yxy*h ) + 
                      k.xxx*map( p + k.xxx*h ) );
}

 vec4 mainImage(  vec2 fragCoord )
{

    vec2 uv = fragCoord/u_resolution.xy;
    
    // screen size is 6m x 6m
    float screenSize=6.0;
    vec3 rayOri = vec3((uv - 0.5) * vec2(u_resolution.x/u_resolution.y, 1.0) * screenSize, 2.0);
    vec3 rayDir = vec3(0.0, 0.0, -1.0);
    
    float depth = 0.0;
    vec3 p;
    
    for(int i = 0; i < 64; i++) {
        p = rayOri + rayDir * depth;
        float dist = map(p);
        depth += dist;
        if (dist < 1e-7) {
            break;
        }
    }
    
    // background 
    depth = min(4.0, depth);
    vec3 n = calcNormal(p);
    float topBrightness=0.577;
    float b = max(0.0, dot(n, vec3(topBrightness)));
    float colTime = 5.8;
    float colShift=(sin(u_time/10.0)+1.0)/2.0*10.0;
    float saturation = 0.8;
    float whiteBalance = 0.5;
    vec3 colorBase = vec3(0.5,3,4);
    vec3 col = (whiteBalance + saturation * cos((b + colTime+colShift) + uv.xyx * 2.0 + colorBase)) * (0.75 + b * 0.35);
    col *= exp( -depth * 0.15 );
    
    // maximum thickness is 2m in alpha channel
    // return vec4(col, 3.0 - (depth - 0.5) / 2.0);
    return vec4(col, 1.0 );
}

void main() {

    gl_FragColor = mainImage(gl_FragCoord.xy);//vec4(1.0,0.0,1.0,1.0);
}