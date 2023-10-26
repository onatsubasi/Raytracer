#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <limits>


using namespace std;
using namespace parser;
typedef unsigned char RGB[3];

class Ray
{
public:
    Vec3f origin{};
    Vec3f direction{};
    bool isShadow;
    int depth;

    Ray(Vec3f origin, Vec3f direction, bool isShadow = false, int depth = 0){
        this->origin = origin;
        this->direction = direction;
        this-> isShadow = isShadow;
        this->depth = depth;
    }
};

class Hit
{
public:
    bool hit; //returns true if ray hits an object.
    float t;
    Vec3f x; //Intersection point.
    Vec3f normal; // surface normal (normalized)
    int materialID;
    Vec3i color; // color of the intersection point
};

Hit MeshIntersectionCheck(Ray ray, Mesh &mesh, const Scene &scene);

Vec3f CrossProduct(Vec3f first, Vec3f second){
    Vec3f result = {0,0,0};
    result.x = first.y * second.z - first.z * second.y;
    result.y = first.z * second.x - first.x * second.z;
    result.z = first.x * second.y - first.y * second.x;

    return result;
}
Vec3i Clamp(Vec3i &vector)
{
    vector.x = fmin(vector.x, 255);
    vector.y = fmin(vector.y, 255);
    vector.z = fmin(vector.z, 255);
    return vector;
}


float DotProduct(Vec3f first, Vec3f second){
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

Vec3f NegateVector(Vec3f vector)
{
    return {-vector.x, -vector.y, -vector.z};
}

Vec3f SumVectors(Vec3f first, Vec3f second)
{
    return {first.x + second.x, first.y + second.y, first.z + second.z };
}

Vec3f SubtractVectors(Vec3f first, Vec3f second)
{
    return {first.x - second.x, first.y - second.y, first.z - second.z };
}

Vec3f MultiplyVector(Vec3f v, float c)
{
    return {v.x * c, v.y * c, v.z * c};
}

Vec3f Normalize(Vec3f vector)
{
    float magnitude = sqrtf(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
    return {vector.x / magnitude, vector.y / magnitude, vector.z / magnitude };
}

float Determinant(Vec3f v1, Vec3f v2, Vec3f v3)
{
    return v1.x * (v2.y * v3.z - v3.y * v2.z) - v1.y * (v2.x * v3.z - v3.x * v2.z) + v1.z * (v2.x * v3.y - v3.x * v2.y);
}

Vec3f FindNormal(Vec3f a, Vec3f b, Vec3f c)
{
    return Normalize(CrossProduct(SubtractVectors(b,a), SubtractVectors(c,a)));  // May not be true. The direction may be different. need to check.
}

float CalculateDistance(Vec3f a, Vec3f b)
{
    return sqrtf(pow((a.x-b.x),2)+pow((a.y-b.y),2)+pow((a.z-b.z),2));
}

float PixelHeight(const Camera &camera)
{
    // x-> left y-> right z-> bottom w-> top
    return (camera.near_plane.w - camera.near_plane.z)/ camera.image_height;
}

float PixelWidth(const Camera &camera)
{
    return (camera.near_plane.y - camera.near_plane.x)/ camera.image_width;
}

Vec3i FindBackgroundColor(const Scene &scene)
{
    return scene.background_color;
}

Vec3f* BuildPixelPositionArray(const Camera &camera)
{
    float pw = PixelWidth(camera);
    float ph = PixelHeight(camera);

    Vec3f u = CrossProduct(camera.up, NegateVector(camera.gaze));

    // -----Initialize the matrix ----------
    Vec3f* result = new Vec3f[camera.image_width * camera.image_height];

    for (int i =0; i< camera.image_width * camera.image_height; i++)
    {
        result[i] = {0,0,0};
    }

    int x = 0;
    int y = 0;
    for(int i = 0; i < camera.image_width * camera.image_height; i++)
    {
        
        if(y >= camera.image_height)
        {
            y = 0;
            x++;
        }

        // do calculation
        result[i].y = -(camera.position.x +
                                u.x * (camera.near_plane.x + ( x + 0.5f) * pw )+
                                camera.up.x * (camera.near_plane.w - (y + 0.5f)*ph)+
                                camera.gaze.x * camera.near_distance);

        result[i].x = -(camera.position.y +
                                u.y * (camera.near_plane.x + ( x + 0.5f) * pw )+
                                camera.up.y * (camera.near_plane.w - (y + 0.5f)*ph)+
                                camera.gaze.y * camera.near_distance);

        result[i].z = camera.position.z +
                                u.z * (camera.near_plane.x + ( x + 0.5f) * pw )+
                                camera.up.z * (camera.near_plane.w - (y + 0.5f) * ph)+
                                camera.gaze.z * camera.near_distance;

        y++;

    }

    return result;
}

Ray SendRay(const Camera &camera, int i, int j)
{
    // find l,r,t,b
    float left = camera.near_plane.x;
    float right = camera.near_plane.y;
    float bottom = camera.near_plane.z;
    float top = camera.near_plane.w;

    //find pixel positions on the image plane
    float su = (right - left)*(j + 0.5)/camera.image_width;
    float sv = (top - bottom)*(i + 0.5)/camera.image_height;

    Vec3f e = camera.position;
    Vec3f w = NegateVector(camera.gaze);
    Vec3f v = camera.up;
    Vec3f u = Normalize(CrossProduct(v,w));

    Vec3f m{}, q{};

    m.x = e.x - (w.x * camera.near_distance);
    m.y = e.y - (w.y * camera.near_distance);
    m.z = e.z - (w.z * camera.near_distance);


    q.x = m.x + u.x*left + v.x*top;
    q.y = m.y + u.y*left + v.y*top;
    q.z = m.z + u.z*left + v.z*top;

    Vec3f s,d;
    s.x = q.x + u.x*su - v.x * sv;
    s.y = q.y + u.y*su - v.y * sv;
    s.z = q.z + u.z*su - v.z * sv;

    d = Normalize(SubtractVectors(s,e));
    Ray ray {e, d};
    ray.isShadow = false;

    return ray;
}

Hit SphereIntersection(Ray ray, const Sphere &sphere, Vec3f center)
{
    Hit result_hit{};
    double A,B,C,t;
    Vec3f o_minus_c = SubtractVectors(ray.origin, center);
    A = DotProduct(ray.direction,ray.direction);
    B = 2 * DotProduct(ray.direction, o_minus_c);
    C = DotProduct(o_minus_c, o_minus_c) - (sphere.radius * sphere.radius);

    double delta = (B*B) - (4 * A * C);

    if(delta < 0.0)
    {
        result_hit.hit = false;
        return result_hit;
    }

    if(delta == 0.0f)
    {
        t = (-B)/(2 * A);
        Vec3f intersectionPoint = SumVectors(ray.origin, MultiplyVector(ray.direction, t));
        result_hit.hit = true,
        result_hit.t = t;
        result_hit.materialID = sphere.material_id;
        result_hit.x = intersectionPoint;
        result_hit.normal = SubtractVectors(intersectionPoint, center);
        return result_hit;
    }

    double t1 = (-B + sqrtf(delta))/(2 * A);
    double t2 = (-B - sqrtf(delta))/(2 * A);
    if(t1 < 0 && t2 < 0)
    {
        result_hit.hit = false;
        return result_hit;
    }
    if(t1 > 0 && t2 > 0)    
    {
        t = fmin(t1,t2);
    }
    else
    {
        t = fmax(t1,t2);
    }    

    Vec3f intersectionPoint = SumVectors(ray.origin, MultiplyVector(ray.direction, t));
    result_hit.hit = true;
    result_hit.t = t;
    result_hit.materialID = sphere.material_id;
    result_hit.x = intersectionPoint;
    result_hit.normal = Normalize(SubtractVectors(intersectionPoint, center));

    return result_hit;
}

Hit TriangleIntersection(const Ray &ray, Vec3f a, Vec3f b, Vec3f c)
{
    Hit result_hit{};
    result_hit.hit = false;

    float alpha, beta, gamma, t, detA;
    Vec3f minusD = NegateVector(ray.direction);
    Vec3f a_minus_c = SubtractVectors(a,c);
    Vec3f b_minus_c = SubtractVectors(b,c);
    Vec3f o_minus_c = SubtractVectors(ray.origin,c);
    // ---Calculate Determinant of A---
    detA = Determinant(a_minus_c, b_minus_c, minusD);
    if(detA == 0.0)
    {
        result_hit.hit = false;
        return result_hit;
    }

    alpha = Determinant(o_minus_c, b_minus_c, minusD) / detA;
    beta = Determinant(a_minus_c, o_minus_c, minusD) / detA;
    t = Determinant(a_minus_c, b_minus_c, o_minus_c) / detA;
    gamma = 1.0 - (alpha + beta);

    if((0 <= alpha && alpha <= 1) && (0 <= beta && beta <= 1) && (0 <= gamma && gamma <= 1) && t > 0)
    {
        Vec3f intersectionPoint = SumVectors(ray.origin, MultiplyVector(ray.direction, t));
        result_hit.hit = true;
        result_hit.t = t;
        result_hit.x = intersectionPoint;
        result_hit.normal = FindNormal(a,b,c);
        return result_hit;
    }

    result_hit.hit = false;
    return result_hit;
}

Hit MeshIntersection(Ray ray, const Mesh &mesh, const Scene &scene)
{
    int faceCount = mesh.faces.size();
    float t_minMesh = numeric_limits<float>::infinity();
    Hit resultHit{};

    for(int faceIndex = 0; faceIndex < faceCount; faceIndex++)
    {
        Vec3f a = scene.vertex_data[mesh.faces[faceIndex].v0_id - 1];
        Vec3f b = scene.vertex_data[mesh.faces[faceIndex].v1_id - 1];
        Vec3f c = scene.vertex_data[mesh.faces[faceIndex].v2_id - 1];

        Hit hit = TriangleIntersection(ray, a, b, c);
        if(hit.hit && hit.t < t_minMesh)
        {
            t_minMesh = hit.t;
            hit.materialID = mesh.material_id;
            resultHit = hit;
        }
    }
    return resultHit;
}


Hit MeshIntersectionCheck(Ray ray, const Mesh &mesh, const Scene &scene)
{
    int faceCount = mesh.faces.size();
    float t_minMesh = numeric_limits<float>::infinity();
    Hit resultHit{};

    for(int faceIndex = 0; faceIndex < faceCount; faceIndex++)
    {
        Vec3f a = scene.vertex_data[mesh.faces[faceIndex].v0_id - 1];
        Vec3f b = scene.vertex_data[mesh.faces[faceIndex].v1_id - 1];
        Vec3f c = scene.vertex_data[mesh.faces[faceIndex].v2_id - 1];

        Hit hit = TriangleIntersection(ray, a, b, c);
        if(hit.hit)
        {
            t_minMesh = hit.t;
            hit.materialID = mesh.material_id;
            resultHit = hit;
            return resultHit;
        }
    }
    return resultHit;
}

Vec3f CalculateAmbientShading(const Scene &scene, Material mat)
{
    Vec3f newColor{};

    newColor.x = scene.ambient_light.x * mat.ambient.x;
    newColor.y = scene.ambient_light.y * mat.ambient.y;
    newColor.z = scene.ambient_light.z * mat.ambient.z;

    return newColor;
}

Vec3i ComputeColor(Hit hit, Ray ray,Scene scene);

Vec3i CalculateShading(const Scene &scene, Material mat, Hit hit, Ray viewRay)
{
    int numOfSpheres = scene.spheres.size();
    int numOfTriangles = scene.triangles.size();
    int numOfMeshes = scene.meshes.size();
    int numOfLights = scene.point_lights.size();
    Vec3i newColor = {0,0,0};

    if(mat.is_mirror)
    {
        /***
        Vec3f reflectionDirection = MultiplyVector((SubtractVectors(hit.normal,viewRay.direction)),2*(DotProduct(hit.normal,viewRay.direction)));

        //Ray reflectionRay = Ray(hit.x,reflectionDirection); // reflection ray
        Ray reflectionRay {hit.x, reflectionDirection};
        reflectionRay.depth++;
        Vec3i computedColor = ComputeColor(hit,reflectionRay,scene);

        newColor.x += computedColor.x * mat.mirror.x;
        newColor.y += computedColor.y * mat.mirror.y;
        newColor.z += computedColor.z * mat.mirror.z;
        ***/
        
    }

    for(int light_index = 0; light_index < numOfLights; light_index++)
    {
        PointLight light = scene.point_lights[light_index];
        bool isClear = true;

        Vec3f epsilonedPoint = SumVectors(MultiplyVector(hit.normal, scene.shadow_ray_epsilon), hit.x);
        float dist = CalculateDistance(light.position,hit.x);

        Vec3f l = Normalize(SubtractVectors(light.position,hit.x));

        Ray currentRay = {epsilonedPoint, l};
        Hit currentHit {};
        // float t_min = numeric_limits<float>::infinity();
        // ---Check for intersections
        for(int t_index = 0; t_index < numOfTriangles && isClear; t_index++)
        {
            Face points = scene.triangles[t_index].indices;
            int pointAIndex = points.v0_id - 1;
            int pointBIndex = points.v1_id - 1;
            int pointCIndex = points.v2_id - 1;
            Hit newHit = TriangleIntersection(currentRay,scene.vertex_data[pointAIndex],scene.vertex_data[pointBIndex],scene.vertex_data[pointCIndex]);

            if(newHit.hit && newHit.t < dist)
            {
                isClear = false;
            }
        }

        for(int s_index = 0; s_index < numOfSpheres && isClear; s_index++)
        {
            int centerID = scene.spheres[s_index].center_vertex_id - 1;
            Hit newHit = SphereIntersection(currentRay,scene.spheres[s_index],scene.vertex_data[centerID]);

            if(newHit.hit && newHit.t < dist)
            {
                isClear = false;
            }
        }

        for(int m_index = 0; m_index < numOfMeshes && isClear; m_index++)
        {
            Hit newHit = MeshIntersectionCheck(currentRay, scene.meshes[m_index], scene);

            if(newHit.hit && newHit.t < dist)
            {
                isClear = false;
            }
        }

        if(!isClear)
            continue; //try the other light

        float cosTheta = DotProduct(hit.normal,l);
        if(cosTheta < 0)
        {
            continue;
        }

        Vec3f h = Normalize(SubtractVectors(l,viewRay.direction));
        float cosAlpha = fmax(0,DotProduct(h,hit.normal));

        Vec3f L_i = {light.intensity.x / (dist * dist), light.intensity.y / (dist * dist), light.intensity.z / (dist * dist)};
        //Vec3f color = CalculateAmbientShading(scene,mat);

        newColor.x += L_i.x * mat.diffuse.x * cosTheta;
        newColor.x += L_i.x * mat.specular.x * powf(cosAlpha, mat.phong_exponent);
        //newColor.x += color.x;


        newColor.y += L_i.y * mat.diffuse.y * cosTheta;
        newColor.y += L_i.y * mat.specular.y * powf(cosAlpha, mat.phong_exponent);
        //newColor.y += color.y;

        newColor.z += L_i.z * mat.diffuse.z * cosTheta;
        newColor.z += L_i.z * mat.specular.z * powf(cosAlpha, mat.phong_exponent);
        //newColor.z += color.z;
    }

    // add ambient color just once
    Vec3f ambientColor = CalculateAmbientShading(scene,mat);
    newColor.x += ambientColor.x;
    newColor.y += ambientColor.y;
    newColor.z += ambientColor.z;

    return Clamp(newColor);
}

Vec3i ComputeColor(Hit hit, Ray ray,Scene scene)
{
    if(ray.depth > scene.max_recursion_depth)
    {
        return {0,0,0};
    }

    if(hit.hit)
    {
        return CalculateShading(scene,scene.materials[hit.materialID-1],hit,ray);
    }
    else if(ray.depth == 0)
    {
        return scene.background_color;
    }
    
    return {0,0,0};
}


int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

    const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };

    Vec3i backGroundColor = FindBackgroundColor(scene);
    Camera currentCam = scene.cameras[0];
    int width = currentCam.image_width , height = currentCam.image_height;
   
   int numOfSpheres = scene.spheres.size();
   int numOfTriangles = scene.triangles.size();
   int numOfMeshes = scene.meshes.size();

    unsigned char* image = new unsigned char [width * height * 3];
    Vec3f* pixelPositions = BuildPixelPositionArray(currentCam);

    int i = 0;
    int j = 0; // For tracking current pixel ID
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            //int colIdx = x / columnWidth;
            //image[i++] = BAR_COLOR[colIdx][0];
            //image[i++] = BAR_COLOR[colIdx][1];
            //image[i++] = BAR_COLOR[colIdx][2];

            //Ray currentRay = Ray(currentCam.position, Normalize(SubtractVectors(pixelPositions[j],currentCam.position)));
            Ray currentRay = SendRay(currentCam, y, x);
            Hit currentHit{};
            float t_min = std::numeric_limits<float>::infinity();

            for(int s_index = 0; s_index < numOfSpheres; s_index++)
            {
                int centerID = scene.spheres[s_index].center_vertex_id - 1;
                Hit newHit = SphereIntersection(currentRay,scene.spheres[s_index],scene.vertex_data[centerID]);

                if(newHit.hit && (newHit.t < t_min))
                {
                    t_min = newHit.t;
                    currentHit = newHit;
                    currentHit.materialID = scene.spheres[s_index].material_id;
                }

            }

            for(int t_index = 0; t_index < numOfTriangles; t_index++)
            {
                Face points = scene.triangles[t_index].indices;
                int pointAIndex = points.v0_id - 1;
                int pointBIndex = points.v1_id - 1;
                int pointCIndex = points.v2_id - 1;
                Hit newHit = TriangleIntersection(currentRay,scene.vertex_data[pointAIndex],scene.vertex_data[pointBIndex],scene.vertex_data[pointCIndex]);

                if(newHit.hit && (newHit.t < t_min))
                {
                    t_min = newHit.t;
                    currentHit = newHit;
                    currentHit.materialID = scene.triangles[t_index].material_id;
                }
            }

            for(int m_index = 0; m_index < numOfMeshes; m_index++)
            {
                Hit newHit = MeshIntersection(currentRay, scene.meshes[m_index], scene);

                if(newHit.hit && (newHit.t < t_min))
                {
                    t_min = newHit.t;
                    currentHit = newHit;
                    currentHit.materialID = scene.meshes[m_index].material_id;
                }
            }
//            int centerID = scene.spheres[0].center_vertex_id - 1;
//            Hit currentHit = SphereIntersection(currentRay,scene.spheres[0],scene.vertex_data[centerID]);
            
            if(currentHit.hit)
            {
                Vec3i color = CalculateShading(scene, scene.materials[currentHit.materialID-1],currentHit, currentRay);

                image[i++] =  round(color.x); // R
                image[i++] =  round(color.y); // G
                image[i++] =  round(color.z); // B
            }
            else
            {
                image[i++] = backGroundColor.x;
                image[i++] = backGroundColor.y;
                image[i++] = backGroundColor.z;
            }

            j++; // increment pixelID
        }
    }

    write_ppm("test.ppm", image, width, height);

}
