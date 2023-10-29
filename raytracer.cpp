#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <limits>
#include <chrono>
#include <thread>
#include <mutex>

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

Vec3f CrossProduct(const Vec3f &first, const Vec3f &second){
    Vec3f result = {0,0,0};
    result.x = first.y * second.z - first.z * second.y;
    result.y = first.z * second.x - first.x * second.z;
    result.z = first.x * second.y - first.y * second.x;

    return result;
}
Vec3f Clamp(const Vec3f &vector)
{
    Vec3f result{};
    result.x = fminf(vector.x, 255);
    result.y = fminf(vector.y, 255);
    result.z = fminf(vector.z, 255);
    return result;
}


float DotProduct(const Vec3f &first, const Vec3f &second){
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

Vec3f NegateVector(const Vec3f &vector)
{
    return {-vector.x, -vector.y, -vector.z};
}

Vec3f SumVectors(const Vec3f &first, const Vec3f &second)
{
    return {first.x + second.x, first.y + second.y, first.z + second.z };
}

Vec3f SubtractVectors(const Vec3f &first, const Vec3f &second)
{
    return {first.x - second.x, first.y - second.y, first.z - second.z };
}

Vec3f MultiplyVector(const Vec3f &v, float c)
{
    return {v.x * c, v.y * c, v.z * c};
}

Vec3f Normalize(const Vec3f &vector)
{
    float magnitude = sqrtf(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
    return {vector.x / magnitude, vector.y / magnitude, vector.z / magnitude };
}

float Determinant(const Vec3f &v1, const Vec3f &v2, const Vec3f &v3)
{
    return v1.x * (v2.y * v3.z - v3.y * v2.z) - v1.y * (v2.x * v3.z - v3.x * v2.z) + v1.z * (v2.x * v3.y - v3.x * v2.y);
}

Vec3f FindNormal(const Vec3f &a, const Vec3f &b, const Vec3f &c)
{
    return Normalize(CrossProduct(SubtractVectors(b,a), SubtractVectors(c,a)));
}

float CalculateDistance(const Vec3f &a, const Vec3f &b)
{
    return sqrtf(powf((a.x-b.x),2)+powf((a.y-b.y),2)+powf((a.z-b.z),2));
}

Vec3f** CalculateJaggedMeshNormals(Scene &scene)
{
    int rows = scene.meshes.size();
    Vec3f** resultArray = new Vec3f*[rows];

    for(int i = 0; i < rows; i++)
    {
        int numOfFaces = scene.meshes[i].faces.size();
        resultArray[i] = new Vec3f[numOfFaces];

        for(int j = 0; j < numOfFaces; j++)
        {
            Face currentFace = scene.meshes[i].faces[j];
            Vec3f currentNormal = FindNormal(scene.vertex_data[currentFace.v0_id-1],scene.vertex_data[currentFace.v1_id-1],scene.vertex_data[currentFace.v2_id-1]);
            resultArray[i][j] = currentNormal;
        }
    }

    return resultArray;
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

Hit SphereIntersection(const Ray &ray, const Sphere &sphere, const Vec3f &center)
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

Hit TriangleIntersection(const Ray &ray, const Vec3f &a, const Vec3f &b, const Vec3f &c, bool isMesh, int meshID, int faceID, Vec3f** &faceNormals)
{
    Hit result_hit{};
    result_hit.hit = false;

    Vec3f faceNormal;
    if(!isMesh)
    {
        faceNormal = FindNormal(a,b,c);
    }
    else
    {
        faceNormal = faceNormals[meshID][faceID];
    }
    
    if(!ray.isShadow && DotProduct(faceNormal, ray.direction) > 0)
    {
        return result_hit;
    }
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
        result_hit.normal = faceNormal;
        return result_hit;
    }

    result_hit.hit = false;
    return result_hit;
}

Hit MeshIntersection(const Ray &ray, const Mesh &mesh, const Scene &scene, int meshID, Vec3f** &faceNormals, double dist = 0.0)
{
    int faceCount = mesh.faces.size();
    double t_minMesh = ray.isShadow? dist: numeric_limits<double>::infinity();
    Hit resultHit{};

    for(int faceIndex = 0; faceIndex < faceCount; faceIndex++)
    {
        Vec3f a = scene.vertex_data[mesh.faces[faceIndex].v0_id - 1];
        Vec3f b = scene.vertex_data[mesh.faces[faceIndex].v1_id - 1];
        Vec3f c = scene.vertex_data[mesh.faces[faceIndex].v2_id - 1];

        Hit hit = TriangleIntersection(ray, a, b, c, true, meshID, faceIndex, faceNormals);

        if(hit.hit && hit.t < t_minMesh)
        {
            if(ray.isShadow)
            {
                resultHit = hit;
                return resultHit;
            }
            t_minMesh = hit.t;
            hit.materialID = mesh.material_id;
            resultHit = hit;
        }
    }
    return resultHit;
}

Hit ClosestHit(const Ray &ray, const Scene &scene, Vec3f** &faceNormals)
{
    int numOfSpheres = scene.spheres.size();
    int numOfTriangles = scene.triangles.size();
    int numOfMeshes = scene.meshes.size();
    Hit currentHit{};
    float t_min = std::numeric_limits<float>::infinity();

    for(int s_index = 0; s_index < numOfSpheres; s_index++)
    {
        int centerID = scene.spheres[s_index].center_vertex_id - 1;
        Hit newHit = SphereIntersection(ray,scene.spheres[s_index],scene.vertex_data[centerID]);

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
        Hit newHit = TriangleIntersection(ray,scene.vertex_data[pointAIndex],scene.vertex_data[pointBIndex],scene.vertex_data[pointCIndex], false, 0,0, faceNormals);

        if(newHit.hit && (newHit.t < t_min))
        {
            t_min = newHit.t;
            currentHit = newHit;
            currentHit.materialID = scene.triangles[t_index].material_id;
        }
    }

    for(int m_index = 0; m_index < numOfMeshes; m_index++)
    {
        Hit newHit = MeshIntersection(ray, scene.meshes[m_index], scene, m_index, faceNormals);

        if(newHit.hit && (newHit.t < t_min))
        {
            t_min = newHit.t;
            currentHit = newHit;
            currentHit.materialID = scene.meshes[m_index].material_id;
        }
    }
    return currentHit;
}

bool InShadow(const Hit &hit, const Scene &scene, const PointLight &I, Vec3f** &faceNormals)
{
    int numOfSpheres = scene.spheres.size();
    int numOfTriangles = scene.triangles.size();
    int numOfMeshes = scene.meshes.size();

    Vec3f lightPos = I.position;
    Vec3f epsilonedPoint = SumVectors(MultiplyVector(hit.normal, scene.shadow_ray_epsilon), hit.x);
    double dist = CalculateDistance(lightPos,epsilonedPoint);
    Vec3f direction = Normalize(SubtractVectors(lightPos, epsilonedPoint));
    Ray shadowRay {epsilonedPoint, direction, true, 0}; 

    for(int t_index = 0; t_index < numOfTriangles; t_index++)
    {
        Face points = scene.triangles[t_index].indices;
        int pointAIndex = points.v0_id - 1;
        int pointBIndex = points.v1_id - 1;
        int pointCIndex = points.v2_id - 1;
        Hit newHit = TriangleIntersection(shadowRay, scene.vertex_data[pointAIndex],scene.vertex_data[pointBIndex],scene.vertex_data[pointCIndex], false, 0,0,faceNormals);
        if(newHit.hit && newHit.t < dist) return true;
    }
    for(int s_index = 0; s_index < numOfSpheres; s_index++)
    {
        Sphere sphere = scene.spheres[s_index];
        int centerID = scene.spheres[s_index].center_vertex_id - 1;
        Hit newHit = SphereIntersection(shadowRay,sphere,scene.vertex_data[centerID]);
        if(newHit.hit && newHit.t < dist) return true;
    }

    for(int m_index = 0; m_index < numOfMeshes; m_index++)
    {
        Mesh mesh = scene.meshes[m_index];
        if(MeshIntersection(shadowRay, mesh, scene, m_index, faceNormals, dist).hit)
            return true;
    }
    return false;
}
Vec3f ApplyShading(const Ray &viewRay, const Hit &hit, const Scene &scene, Vec3f** &faceNormals);
Vec3f ComputeColor(const Ray &viewRay, const Scene &scene, const Hit &closestHit, Vec3f** &faceNormals);

Ray Reflect(const Ray &ray, const Hit &hit, const float epsilon)
{
    Vec3f epsilonedPoint = SumVectors(MultiplyVector(hit.normal, epsilon), hit.x);
    Vec3f wo = NegateVector(ray.direction);
    Vec3f wr = Normalize(SubtractVectors(MultiplyVector(hit.normal,(2 * DotProduct(hit.normal, wo))) , wo)); //TODO need to check if something goes wrong.
    return {epsilonedPoint, wr};
}

Vec3f DiffuseTerm(const Hit &hit, const PointLight &light, const Material &material)
{
    Vec3f result{};
    Vec3f l = Normalize(SubtractVectors(light.position, hit.x));
    float cosTheta = DotProduct(hit.normal, l);
    if(cosTheta <= 0) return {0,0,0};
    float dist = CalculateDistance(light.position, hit.x);
    Vec3f L_i = {light.intensity.x / (dist * dist), light.intensity.y / (dist * dist), light.intensity.z / (dist * dist)};
    result.x = L_i.x * material.diffuse.x * cosTheta;
    result.y = L_i.y * material.diffuse.y * cosTheta;
    result.z = L_i.z * material.diffuse.z * cosTheta;
    return result;
}

Vec3f SpecularTerm(const Ray &viewRay, const Hit &hit, const PointLight &light, const Material &material)
{
    Vec3f result{};
    Vec3f l = Normalize(SubtractVectors(light.position, hit.x));
    Vec3f h = Normalize(SubtractVectors(l,viewRay.direction));
    float cosTheta = DotProduct(hit.normal, l);
    if(cosTheta <= 0) return {0,0,0};
    float cosAlpha = fmax(0, DotProduct(h,hit.normal));
    float dist = CalculateDistance(light.position, hit.x);
    Vec3f L_i = {light.intensity.x / (dist * dist), light.intensity.y / (dist * dist), light.intensity.z / (dist * dist)};
    float phong = powf(cosAlpha,material.phong_exponent);
    result.x = L_i.x * material.specular.x * phong;
    result.y = L_i.y * material.specular.y * phong;
    result.z = L_i.z * material.specular.z * phong;
    return result;
}

Vec3f ApplyShading(const Ray &viewRay, const Hit &hit, const Scene &scene, Vec3f** &faceNormals)
{
    int numOfLights = scene.point_lights.size();
    Material material = scene.materials[hit.materialID - 1];
    Vec3f color {material.ambient.x * scene.ambient_light.x, material.ambient.y * scene.ambient_light.y, material.ambient.z * scene.ambient_light.z};
    if(material.is_mirror)
    {
        Ray reflectionRay = Reflect(viewRay, hit, scene.shadow_ray_epsilon);
        reflectionRay.depth = viewRay.depth + 1;
        Hit closestHit = ClosestHit(reflectionRay, scene, faceNormals);
        Vec3f tempColor = ComputeColor(reflectionRay, scene, closestHit, faceNormals);

        color.x += tempColor.x * material.mirror.x;
        color.y += tempColor.y * material.mirror.y;
        color.z += tempColor.z * material.mirror.z;
    }
    for(int l_index = 0; l_index < numOfLights; l_index++)
    {
        PointLight I = scene.point_lights[l_index];
        if(!InShadow(hit, scene, I, faceNormals))
        {

            Vec3f diffuseTerm = DiffuseTerm(hit, I, material);
            Vec3f specularTerm = SpecularTerm(viewRay, hit, I, material);
            color.x += diffuseTerm.x + specularTerm.x;
            color.y += diffuseTerm.y + specularTerm.y;
            color.z += diffuseTerm.z + specularTerm.z;
        }
    }
    return color;
}

Vec3f ComputeColor(const Ray &viewRay, const Scene &scene, const Hit &closestHit, Vec3f** &faceNormals)
{
    if(viewRay.depth > scene.max_recursion_depth)
        return {0,0,0};

    if(closestHit.hit)
        return ApplyShading(viewRay, closestHit, scene, faceNormals);
    else if (viewRay.depth == 0)
    {
        Vec3f result{};
        result.x = scene.background_color.x;
        result.y = scene.background_color.y;
        result.z = scene.background_color.z;
        return result;
    }
    else
        return {0,0,0};
}

void f1(Scene &scene, Camera &currentCam, int width, unsigned char* &image, Vec3f** JaggedArrayMeshNormals, int y, int x)
{
    int pixel_index = (y * width + x) * 3;
    Ray currentRay = SendRay(currentCam, y, x);
    currentRay.depth = 0;
    Hit closestHit = ClosestHit(currentRay, scene, JaggedArrayMeshNormals);
    Vec3f color = Clamp(ComputeColor(currentRay, scene, closestHit, JaggedArrayMeshNormals));

    image[pixel_index] =  round(color.x); // R
    image[pixel_index+1] =  round(color.y); // G
    image[pixel_index+2] =  round(color.z); // B
}

int main(int argc, char* argv[])
{
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    Vec3f** JaggedArrayMeshNormals = CalculateJaggedMeshNormals(scene);

    int numOfCameras = scene.cameras.size();
    for(int c_index = 0; c_index < numOfCameras; c_index++)
    {
        Camera currentCam = scene.cameras[c_index];
        auto start = std::chrono::high_resolution_clock::now();
        cout << "Rendering " << currentCam.image_name.c_str() << " ...\n";
        int width = currentCam.image_width , height = currentCam.image_height;

        unsigned int threadCount = std::thread::hardware_concurrency();
        std::vector<thread> myThreads;
        unsigned char *image = new unsigned char[width * height * 3];

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; x += threadCount) {
                myThreads.clear(); // Clear the vector before launching new threads

                for (int i = 0; i < threadCount && x + i < width; i++)
                {
                    myThreads.emplace_back(f1, std::ref(scene), std::ref(currentCam), width, std::ref(image),
                                               JaggedArrayMeshNormals, y, x + i);
                }


                // Join the threads before moving to the next iteration
                for (std::thread &t : myThreads) {
                    t.join();
                }
            }
        }
        write_ppm(currentCam.image_name.c_str(), image, width, height);
        cout << "Rendering Complete.\n";
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        cout << "Time taken to render " << currentCam.image_name.c_str() << ": " << (int) duration.count()/60 << " min " << fmod(duration.count(),60) << " sec.\n";
    }
}
