// Tampere University
// COMP.CE.430 Computer Graphics Coding Assignment 2 (2021)
//
// Write your name and student id here:
//   Max Gratschew, H283272
//
// Mark here with an X which functionalities you implemented.
// Note that different functionalities are worth different amount of points.
//
// Name of the functionality      |Done| Notes
//-------------------------------------------------------------------------------
// example functionality          | X  | Example note: control this with var YYYY
// Mandatory functionalities ----------------------------------------------------
//   Perspective projection       |  X | 
//   Phong shading                |  X | 
//   Camera movement and rotation |  X | 
//   Sharp shadows                |  X | comment #define SOFT_SHADOWS from the line 48
// Extra functionalities --------------------------------------------------------
//   Tone mapping                 |  X | uncomment #define TONE_MAPPING from the line 51
//   PBR shading                  |    | 
//   Soft shadows                 |  X | 
//   Sharp reflections            |  X | 
//   Glossy reflections           |    | 
//   Refractions                  |    | 
//   Caustics                     |    | 
//   SDF Ambient Occlusions       |  X | 
//   Texturing                    |  X | Put on the sphere as a material
//   Simple game                  |    | 
//   Progressive path tracing     |    | 
//   Basic post-processing        |    | 
//   Advanced post-processing     |    | 
//   Screen space reflections     |    | 
//   Screen space AO              |    | 
//   Simple own SDF               |  X | Torus
//   Advanced own SDF             |  X | Double pyramid
//   Animated SDF                 |  X | Double pyramid's morphing animation
//   Other?                       |    | 


#ifdef GL_ES
precision mediump float;
#endif

#define PI 3.14159265359
#define EPSILON 0.00001
    
// UNCOMMENT THIS TO TOGGLE OFF SOFT_SHADOWS (ENABLING HARD SHADOWS)    
#define SOFT_SHADOWS

// UNCOMMENT THIS TO TOGGLE ON TONE MAPPING
//#define TONE_MAPPING

// LIGHTING
const vec3 AMBIENT = vec3(0.202,0.203,0.205);
const vec3 DIFFUSE = vec3(0.743,0.800,0.792);
const vec3 SPECULAR = vec3(0.670,0.632,0.613);
const float DIFFUSE_INTENSITY = 0.79;
const float SPECULAR_INTENSITY = 0.70;
const float SHININESS = 30.00;

// These definitions are tweakable.

/* Minimum distance a ray must travel. Raising this value yields some performance
 * benefits for secondary rays at the cost of weird artefacts around object
 * edges.
 */
#define MIN_DIST 0.08
/* Maximum distance a ray can travel. Changing it has little to no performance
 * benefit for indoor scenes, but useful when there is nothing for the ray
 * to intersect with (such as the sky in outdoors scenes).
 */
#define MAX_DIST 20.0
/* Maximum number of steps the ray can march. High values make the image more
 * correct around object edges at the cost of performance, lower values cause
 * weird black hole-ish bending artefacts but is faster.
 */
#define MARCH_MAX_STEPS 128
/* Typically, this doesn't have to be changed. Lower values cause worse
 * performance, but make the tracing stabler around slightly incorrect distance
 * functions.
 * The current value merely helps with rounding errors.
 */
#define STEP_RATIO 0.97
/* Determines what distance is considered close enough to count as an
 * intersection. Lower values are more correct but require more steps to reach
 * the surface
 */
#define HIT_RATIO 0.001

// Resolution of the screen
uniform vec2 u_resolution;

// Mouse coordinates
uniform vec2 u_mouse;

// Time since startup, in seconds
uniform float u_time;

struct material
{
    // The color of the surface
    vec4 color;
    // You can add your own material features here!
    
    float type;
    
};

vec2 random(vec2 st){
  st=vec2(dot(st,vec2(127.1,311.7)),dot(st,vec2(269.5,183.3)));
  return-1.+2.*fract(sin(st)*43758.5453123);
}
float noise(vec2 st){
  vec2 i=floor(st);
  vec2 f=fract(st);
  
  vec2 u=f*f*(3.6 -2.6*f);//*(3.-2.*f);
  
  return mix(mix(dot(random(i+vec2(0.,0.)),f-vec2(0.,0.)),
  	dot(random(i+vec2(1.,0.)),f-vec2(1.,0.)),u.x),
  	mix(dot(random(i+vec2(0.,1.)),f-vec2(0.,1.)),
  	dot(random(i+vec2(1.,1.)),f-vec2(1.,1.)),u.x),u.y);
}



// Good resource for finding more building blocks for distance functions:
// https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm

/* Basic box distance field.
 *
 * Parameters:
 *  p   Point for which to evaluate the distance field
 *  b   "Radius" of the box
 *
 * Returns:
 *  Distance to the box from point p.
 */
float box(vec3 p, vec3 b)
{
    vec3 d = abs(p) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

// Rotation matrix around the X axis.
mat3 rotateX(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(1, 0, 0),
        vec3(0, c, -s),
        vec3(0, s, c)
    );
}

// Rotation matrix around the Y axis.
mat3 rotateY(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, 0, s),
        vec3(0, 1, 0),
        vec3(-s, 0, c)
    );
}

/* Rotates point around origin along the X axis.
 *
 * Parameters:
 *  p   The point to rotate
 *  a   The angle in radians
 *
 * Returns:
 *  The rotated point.
 */
vec3 rot_x(vec3 p, float a)
{
    float s = sin(a);
    float c = cos(a);
    return vec3(
        p.x,
        c*p.y-s*p.z,
        s*p.y+c*p.z
    );
}

/* Rotates point around origin along the Y axis.
 *
 * Parameters:
 *  p   The point to rotate
 *  a   The angle in radians
 *
 * Returns:
 *  The rotated point.
 */
vec3 rot_y(vec3 p, float a)
{
    float s = sin(a);
    float c = cos(a);
    return vec3(
        c*p.x+s*p.z,
        p.y,
        -s*p.x+c*p.z
    );
}

/* Rotates point around origin along the Z axis.
 *
 * Parameters:
 *  p   The point to rotate
 *  a   The angle in radians
 *
 * Returns:
 *  The rotated point.
 */
vec3 rot_z(vec3 p, float a)
{
    float s = sin(a);
    float c = cos(a);
    return vec3(
        c*p.x-s*p.y,
        s*p.x+c*p.y,
        p.z
    );
}

/* Each object has a distance function and a material function. The distance
 * function evaluates the distance field of the object at a given point, and
 * the material function determines the surface material at a point.
 */

float sdPyramid_distance( in vec3 p, in float h )
{
    float m2 = h*h + 1.25;
    p.x = p.x-2.0;
    p.y = p.y+0.0;
    p.z = p.z-3.50;
  
    // double pyramid
    p=abs(p);
    
    p.xz = (p.z>p.x) ? p.zx : p.xz;
    p.xz -= 0.1;
	
    // project into face plane (2D)
    vec3 q = vec3( p.z, h*p.y-0.5*p.x, h*p.x+0.5*p.y);
        
    float s = max(-q.x,0.0);
    float t = clamp( (q.y-0.5*q.x)/(m2+0.25), 0.0, 1.0 );
    
    float a = m2*(q.x+s)*(q.x+s) + q.y*q.y;
	float b = m2*(q.x+0.5*t)*(q.x+0.5*t) + (q.y-m2*t)*(q.y-m2*t);
    
    float d2 = max(-q.y,q.x*m2+q.y*0.5) < 0.0 ? 0.0 : min(a,b);
    
    // recover 3D and scale,add sign
    return sqrt( (d2+q.z*q.z)/m2 ) * sign(max(q.z,-p.y));;
}

material pyramid_material( vec3 p ){
    material mat;
    mat.color = vec4(0.0, 2.0, 1., 0.0);
    return mat;
}
float sdTorus_distance( vec3 p )
{
    float le = .2;
    float r1 = .2;
    float r2 = .15;
  	vec3 q = vec3( p.x-1.0, max(abs(p.y+2.5)-le,.00), p.z-2.0 );
  	return length(vec2(length(q.xy)-r1,q.z)) - r2;
}

material sdTorus_material( vec3 p )
{
    material mat;
    mat.color = vec4(1.0, 0.5, 0.9, 0.0);

    
    return mat;
}

float blob_distance(vec3 p)
{
    vec3 q = p - vec3(-0.5, -2.2 + abs(sin(u_time*3.0)), 2.0);
    return length(q) - 0.8 + sin(10.0*q.x)*sin(10.0*q.y)*sin(10.0*q.z)*0.07;
}

material blob_material(vec3 p)
{
    material mat;
    mat.color = vec4(1.0, 0.5, 0.3, 0.0);
    return mat;
}

float sphere_distance(vec3 p)
{
    return length(p - vec3(1.5, -1.8, 4.0)) - 1.2;
}

material sphere_material(vec3 p)
{
    material mat;
    vec2 st = vec2(p.x, p.y);
    
  	st+=noise(st*5.)*abs(2.-sin(u_time*.9))*5.+2.;
  	float splatterFreq=0.7;
  	float redSplatter=smoothstep(.01,.16,noise(st*splatterFreq+10.));
  	float greenSplatter=smoothstep(.01,.12,noise(st*splatterFreq*2.+20.));
  	float blueSplatter=smoothstep(.08,.14,noise(st*splatterFreq*3.+30.));
  
  	vec3 color = vec3(redSplatter, greenSplatter, blueSplatter);
    mat.color = vec4(color, 1.0);
    
    
    return mat;
}

float room_distance(vec3 p)
{
    return max(
        -box(p-vec3(0.0,3.0,3.0), vec3(0.5, 0.5, 0.5)),
        -box(p-vec3(0.0,0.0,0.0), vec3(3.0, 3.0, 6.0))
    );
}

material room_material(vec3 p)
{
    material mat;
    mat.color = vec4(0.945,0.945,0.945,1.000);
    if(p.x <= -2.98) {
        mat.color.rgb = vec3(1.0, 0.0, 0.0);
 
    
    return mat;
    }
    else if(p.x >= 2.98) mat.color.rgb = vec3(0.0, 1.0, 0.0);
    else if( p.x < 3.0 && p.x >-3.0 && p.y < 3.1 && p.y > -3.1 && p.z > 5.9) mat.color.rgb = vec3(1.00,1.00,1.00);
    
    return mat;
}

float crate_distance(vec3 p)
{
    return box(rot_y(p-vec3(-1,-1,5), u_time), vec3(1, 2, 1));
}

material crate_material(vec3 p)
{
    material mat;
    mat.color = vec4(1.0, 1.0, 1.0, 1.0);

    vec3 q = rot_y(p-vec3(-1,-1,5), u_time) * 0.98;
    if(fract(q.x + floor(q.y*2.0) * 0.5 + floor(q.z*2.0) * 0.5) < 0.5)
    {
        mat.color.rgb = vec3(0.0, 1.0, 1.0);
    }
    return mat;
}

/* The distance function collecting all others.
 *
 * Parameters:
 *  p   The point for which to find the nearest surface
 *  mat The material of the nearest surface
 *
 * Returns:
 *  The distance to the nearest surface.
 */
float map(
    in vec3 p,
    out material mat
){
    float min_dist = MAX_DIST*2.0;
    float dist = 0.0;

    dist = blob_distance(p);
    if(dist < min_dist) {
        mat = blob_material(p);
        min_dist = dist;
    }

    dist = room_distance(p);
    if(dist < min_dist) {
        mat = room_material(p);
        min_dist = dist;
    }

    dist = crate_distance(p);
    if(dist < min_dist) {
        mat = crate_material(p);
        min_dist = dist;
    }

    dist = sphere_distance(p);
    if(dist < min_dist) {
        
        // moving splatter texture
        vec2 st = vec2(0.5);
  		st+=noise(st*5.)*abs(2.-sin(u_time*.9))*5.+2.;
  		float splatterFreq=1.;
  		float redSplatter=smoothstep(.01,.16,noise(st*splatterFreq+10.));
  		float greenSplatter=smoothstep(.01,.12,noise(st*splatterFreq*2.+20.));
  		float blueSplatter=smoothstep(.08,.14,noise(st*splatterFreq*3.+30.));
  
  		vec3 color = vec3(redSplatter, greenSplatter, blueSplatter);
        mat = sphere_material(p);
        min_dist = dist;
    }
    
	// own objects here!
    // torus
    dist = sdTorus_distance(p);
    if(dist < min_dist){
        mat = sdTorus_material(p);
        min_dist = dist;
    }
    
    // pyramid
    float rad = 0.2*(0.5-0.45*cos(u_time*3.0));
    float height = 1.5*(0.5+0.45*sin(u_time*1.0));
    dist = sdPyramid_distance(p, height)-rad;
    if(dist < min_dist){
        mat = pyramid_material(p);
        min_dist = dist;
    }
    

    return min_dist;
}

/* Calculates the normal of the surface closest to point p.
 *
 * Parameters:
 *  p   The point where the normal should be calculated
 *  mat The material information, produced as a byproduct
 *
 * Returns:
 *  The normal of the surface.
 *
 * See https://www.iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
 * if you're interested in how this works.
 */
vec3 normal(vec3 p, out material mat)
{
    const vec2 k = vec2(1.0, -1.0);
    return normalize(
        k.xyy * map(p + k.xyy * EPSILON, mat) +
        k.yyx * map(p + k.yyx * EPSILON, mat) +
        k.yxy * map(p + k.yxy * EPSILON, mat) +
        k.xxx * map(p + k.xxx * EPSILON, mat)
    );
}

/* Finds the closest intersection of the ray with the scene.
 *
 * Parameters:
 *  o           Origin of the ray
 *  v           Direction of the ray
 *  max_dist    Maximum distance the ray can travel. Usually MAX_DIST.
 *  p           Location of the intersection
 *  n           Normal of the surface at the intersection point
 *  mat         Material of the intersected surface
 *  inside      Whether we are marching inside an object or not. Useful for
 *              refractions.
 *
 * Returns:
 *  true if a surface was hit, false otherwise.
 */
bool intersect(
    in vec3 o,
    in vec3 v,
    in float max_dist,
    out vec3 p,
    out vec3 n,
    out material mat,
    bool inside
) {
    float t = MIN_DIST;
    float dir = inside ? -1.0 : 1.0;
    bool hit = false;

    for(int i = 0; i < MARCH_MAX_STEPS; ++i)
    {
        p = o + t * v;
        float dist = dir * map(p, mat);
        
        hit = abs(dist) < HIT_RATIO * t;

        if(hit || t > max_dist) break;

        t += dist * STEP_RATIO;
    }

    n = normal(p, mat);

    return hit;
}

//n = surface normal
//rd = view direction
//ld = light direction
vec3 shade(vec3 n, vec3 rd, vec3 ld){
    
    n = normalize(n);
    ld = normalize(ld);
    rd = -rd;
    float lambertian = max(dot(ld, n),0.0);
    float spec = 0.00;
    
    if (lambertian > 0.00) {
        vec3 halfDir = normalize(ld+rd);
        float specAngle = max(dot(halfDir, n),0.0);
        spec = pow(specAngle, SHININESS);
    }
    vec3 color = AMBIENT + lambertian * DIFFUSE*DIFFUSE_INTENSITY + spec* SPECULAR * SPECULAR_INTENSITY;
    return color;
}



float softShadow(vec3 p, vec3 lamp_pos){
	vec3 rd = normalize(lamp_pos -p);
	vec3 ro = p;
    float res = 1.0;
    float tmin = 0.35;
    float tmax = 100.0;

    float t = 1.0;
    material temp_mat;
    for( int i=0; i<16; i++ )
    {
        float h = map( ro + rd*t, temp_mat );
        res = min( res, 1.750*h/t );
        t += clamp( h, 0.02, 0.10 );
        if( h<0.001 || t>tmax ) break;
    }
    return clamp( res, 0.1, 1.0 );
}

float ambientOcclusion(vec3 p, vec3 n) {
    float step = 0.5;
    float ao = 0.0;
    float dist;
    float d;
    material temp_mat;
    for (int i = 1; i <= 3; i++) {
        
        dist = step * float(i);
        d = map(p + n * dist, temp_mat);
  ao += max(0., (dist - d) / dist);  
    }
    
    return 1.0 - ao * 0.1;
}
    
    

/* Calculates the color of the pixel, based on view ray origin and direction.
 *
 * Parameters:
 *  o   Origin of the view ray
 *  v   Direction of the view ray
 *
 * Returns:
 *  Color of the pixel.
 */
vec3 render(vec3 o, vec3 v)
{
    // This lamp is positioned at the hole in the roof.
    vec3 lamp_pos = vec3(0.0, 3.1, 3.0);

    vec3 p, n;
    bool shadow;
    material mat;

    // Compute intersection point along the view ray.
    bool hit = intersect(o, v, MAX_DIST, p, n, mat, false);
    
    // light direction
    vec3 ld = normalize(lamp_pos - p);
    vec3 original_pixel = p;
    vec3 original_normal = n;
    
    // sharp reflection
    vec3 p_r, n_r;
    material reflection_material;
    vec3 reflected_ray = reflect(v,n);
    bool reflection_hit = intersect(p, reflected_ray, MAX_DIST, p_r, n_r, reflection_material, false);
    material reflective_material;
   
    
    map(p, reflective_material);

	if(reflection_hit){
        // white colored materials to work as a mirror
 		if(reflective_material.color.rgb == vec3(1.00,1.00,1.00)){
            mat.color *= reflection_material.color;
            
            // lighting the reflection
        	material temp_mat_r = reflection_material;
        	vec3 n_r2 = normal(p_r, temp_mat_r);
            vec3 ld_r = normalize(lamp_pos - p_r);
        	vec3 lighting = shade(n_r2, v, ld_r);
            
            
            material temp_material;
            vec3 sro = p_r+2.0*EPSILON*n_r;
            float lamp_dist = distance(lamp_pos, sro);
            bool is_in_shadow = intersect(sro, ld_r, lamp_dist, p_r, n_r, temp_material, false);
            if(is_in_shadow){
                
            #ifdef SOFT_SHADOWS
                
			float soft = softShadow(p_r, lamp_pos);
        	mat.color.rgb *= soft;
                
			#else
                
			mat.color.rgb *= 0.3;
                
			#endif
             mat.color.rgb *= AMBIENT;
  
            }
            else{
            mat.color.rgb *= lighting;
            }
            
        }
    }
    
    if(hit){
        
        
        material temp_mat = mat;
        vec3 n = normal(p, temp_mat);
        
      	// phong lighting for later use
        vec3 lighting = shade(n, v, ld);
        
        vec3 sro = p + 2.0 * EPSILON * n;
            
		
        // Shadows
        p = vec3 (0.0);
        float lamp_dist = distance(lamp_pos, sro);
        bool hit = intersect(sro, ld, lamp_dist, p, n, temp_mat, false);
       
        if(hit){
                        
        #ifdef SOFT_SHADOWS
                
		float soft = softShadow(p, lamp_pos);
        mat.color.rgb *= soft;
          
		#else
                
		mat.color.rgb *= 0.3;
                
		#endif
    	mat.color.rgb *= AMBIENT;
            
        }
        
        else {
        mat.color.rgb *= lighting;
        }

        // ambient occlusion
    	float ao = ambientOcclusion(original_pixel, original_normal);
    	mat.color *= ao; 
        
    }

    return mat.color.rgb;
}

/* Calculates the color of the pixel, based on view ray origin and direction.
 *
 * Parameters:
 *  color   material's color
 *
 * Returns:
 *  tone mapped color via simple Reinhard tone mapping
 */
vec3 simpleReinhardToneMapping(vec3 color)
{
    float gamma = 2.2;
	float exposure = 1.5;
	color *= exposure/(1. + color / exposure);
	color = pow(color, vec3(1. / gamma));
	return color;
}

    

void main()
{
    // This is the position of the pixel
    vec2 uv = (gl_FragCoord.xy-.5*u_resolution.xy)/u_resolution.y;
    

    
    float width = u_resolution.x;
    float height = u_resolution.y;

    

    vec2 mouse_pos = u_mouse.xy / u_resolution.xy - 0.5;
    
    // ray origin
    vec3 o = vec3(0.0,0.0, -5.0);
    
 
    // ray direction
    vec3 v = normalize(vec3(uv, 1));
    
    // rotate camera with the mouse
	v *= rotateY(mouse_pos.x) * rotateX(mouse_pos.y);

	#ifdef TONE_MAPPING
    gl_FragColor = vec4(simpleReinhardToneMapping(render(o, v)), 1);
    #else
    gl_FragColor = vec4(render(o, v), 1);
    #endif

}
