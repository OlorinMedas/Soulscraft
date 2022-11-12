#include "fbx_common.h"

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>

#include <unordered_set>

//#define DEBUG_INFO

using namespace std;

struct Point
{
    int x;
    int y;
    int z;

    bool operator==(const Point& other) const
    {
        return (x == other.x
                && y == other.y
                && z == other.z);
    }
};

template <>
struct hash<Point>
{
    size_t operator()(const Point& k) const
    {
        // Compute individual hash values for first, second and third
        // http://stackoverflow.com/a/1646913/126995
        size_t res = 17;
        res = res * 31 + hash<int>()(k.x);
        res = res * 31 + hash<int>()(k.y);
        res = res * 31 + hash<int>()(k.z);
        return res;
    }
};

unordered_set<Point> g_point_cloud;

void segmentPolygon(FbxVector4* p_in, int vertex_count_in,
                    FbxVector4* p_out1, int& vertex_count_out1,
                    FbxVector4* p_out2, int& vertex_count_out2,
                    double pivot, int axis)
{
    float d[12];
    for (int i = 0; i < vertex_count_in; ++i)
    {
        d[i] = pivot - p_in[i][axis];
    }

    int m = 0, n = 0;
    for (int i = 0, j = vertex_count_in - 1; i < vertex_count_in; j = i, ++i)
    {
        bool ina = d[j] >= 0;
        bool inb = d[i] >= 0;
        if (ina != inb)
        {
            float s = d[j] / (d[j] - d[i]);
            p_out1[m][0] = p_in[j][0] + (p_in[i][0] - p_in[j][0]) * s;
            p_out1[m][1] = p_in[j][1] + (p_in[i][1] - p_in[j][1]) * s;
            p_out1[m][2] = p_in[j][2] + (p_in[i][2] - p_in[j][2]) * s;
            p_out2[n] = p_out1[m];
            m++;
            n++;
            // add the i'th point to the right polygon. Do NOT add points that are on the dividing line
            // since these were already added above
            if (d[i] > 0)
            {
                p_out1[m] = p_in[i];
                m++;
            }
            else if (d[i] < 0)
            {
                p_out2[n] = p_in[i];
                n++;
            }
        }
        else // same side
        {
            // add the i'th point to the right polygon. Addition is done even for points on the dividing line
            if (d[i] >= 0)
            {
                p_out1[m] = p_in[i];
                m++;
            }
            if (d[i] <= 0)
            {
                p_out2[n] = p_in[i];
                n++;
            }
            
        }
    }

    vertex_count_out1 = m;
    vertex_count_out2 = n;
}

void segmentPoint(FbxVector4* v0, FbxVector4* v1, double pivot, int axis, FbxVector4& out_segment_point)
{
    double ratio = (pivot - (*v0)[axis]) / ((*v1)[axis] - (*v0)[axis]);
    out_segment_point[axis] = pivot;
    const int second_axis = (axis + 1) % 3;
    out_segment_point[second_axis] = ratio * ((*v1)[second_axis] - (*v0)[second_axis]) + (*v0)[second_axis];
    const int third_axis = (axis + 2) % 3;
    out_segment_point[third_axis] = ratio * ((*v1)[third_axis] - (*v0)[third_axis]) + (*v0)[third_axis];;
}

unsigned int rasterizeTriangle(FbxVector4& v0, FbxVector4& v1, FbxVector4& v2,
                       FbxVector4& lower_bound,
                       FbxVector4& sample_size, FbxVector4& inversed_sample_size,
                       ofstream& output_file)
{
    unsigned int sample_count = 0;
    //calculate triangle bound
    const double max_double = std::numeric_limits<double>::max();
    double x_lower_bound = max_double;
    double x_upper_bound = -max_double;

    FbxVector4 buffer[4][7];
    int vertex_counts[4];
    FbxVector4* vertices, * vertices_row, * p1, * p2, * temp;
    vertices = buffer[0];
    vertices_row = buffer[1];
    p1 = buffer[2];
    p2 = buffer[3];

    vertices[0] = v0;
    vertices[1] = v1;
    vertices[2] = v2;
    vertex_counts[0] = 3;

    for (int i = 0; i < 3; ++i)
    {
        x_lower_bound = FbxMin(x_lower_bound, (vertices[i])[0]);
        x_upper_bound = FbxMax(x_upper_bound, (vertices[i])[0]);
    }

    int x_index_bound[2];
    x_index_bound[0] = (int)::FbxFloor((x_lower_bound - lower_bound[0]) * inversed_sample_size[0]);
    x_index_bound[1] = (int)::FbxFloor((x_upper_bound - lower_bound[0]) * inversed_sample_size[0]) + 1;

    for (int x_index = x_index_bound[0]; x_index < x_index_bound[1]; ++x_index)
    {
        double x_pivot = (x_index + 1) * sample_size[0] + lower_bound[0];
        segmentPolygon(vertices, vertex_counts[0],
                       vertices_row, vertex_counts[1],
                       p1, vertex_counts[2],
                       x_pivot, 0);

        //memcpy(vertices, p1, sizeof(FbxVector4) * 7);
        temp = vertices;
        vertices = p1;
        p1 = temp;
        vertex_counts[0] = vertex_counts[2];

        double z_lower_bound = max_double;
        double z_upper_bound = -max_double;
        for (int i = 0; i < vertex_counts[1]; ++i)
        {
            z_lower_bound = FbxMin(z_lower_bound, vertices_row[i][2]);
            z_upper_bound = FbxMax(z_upper_bound, vertices_row[i][2]);
        }

        int z_index_bound[2];
        z_index_bound[0] = (int)::FbxFloor((z_lower_bound - lower_bound[2]) * inversed_sample_size[2]);
        z_index_bound[1] = (int)::FbxFloor((z_upper_bound - lower_bound[2]) * inversed_sample_size[2]) + 1;

        for (int z_index = z_index_bound[0]; z_index < z_index_bound[1]; ++z_index)
        {
            double z_pivot = (z_index + 1) * sample_size[2] + lower_bound[2];

            segmentPolygon(vertices_row, vertex_counts[1],
                           p1, vertex_counts[2],
                           p2, vertex_counts[3],
                           z_pivot, 2);

            temp = vertices_row;
            vertices_row = p2;
            p2 = temp;
            vertex_counts[1] = vertex_counts[3];

            double y_lower_bound = max_double;
            double y_upper_bound = -max_double;
            for (int i = 0; i < vertex_counts[2]; ++i)
            {
                y_lower_bound = FbxMin(y_lower_bound, p1[i][1]);
                y_upper_bound = FbxMax(y_upper_bound, p1[i][1]);
            }

            int y_index_bound[2];
            y_index_bound[0] = (int)::FbxFloor((y_lower_bound - lower_bound[1]) * inversed_sample_size[1]);
            y_index_bound[1] = (int)::FbxFloor((y_upper_bound - lower_bound[1]) * inversed_sample_size[1]) + 1;

            for (int y_index = y_index_bound[0]; y_index < y_index_bound[1]; ++y_index)
            {
                Point point{x_index, y_index, z_index};
                if (g_point_cloud.count(point) > 0) continue;

                g_point_cloud.insert(point);
                
                output_file.write((char*)&point.x, sizeof(int));
                output_file.write((char*)&point.y, sizeof(int));
                output_file.write((char*)&point.z, sizeof(int));
                sample_count++;
            }
        }
    }

#ifdef DEBUG_INFO
    cout << sample_count << " samples" << endl;
#endif

    return sample_count;
}

void calculateBounds(FbxScene* scene, FbxVector4* bounds)
{
    FbxNode* root_node = scene->GetRootNode();
    int child_count = root_node->GetChildCount();
    for (int index = 0; index < child_count; ++index)
    {
        FbxNode* child_node = root_node->GetChild(index);
        if (child_node == nullptr) continue;
        FbxNodeAttribute* attribute = child_node->GetNodeAttribute();
        FbxNodeAttribute::EType attributeType = attribute->GetAttributeType();
        if (attributeType != FbxNodeAttribute::eMesh) continue;
        FbxDouble3 translation = child_node->LclTranslation.Get();
        FbxDouble3 rotation = child_node->LclRotation.Get();
        FbxDouble3 scaling = child_node->LclScaling.Get();

        translation[0] /= 100.0;
        translation[1] /= 100.0;
        translation[2] /= 100.0;

        scaling[0] /= 100.0;
        scaling[1] /= 100.0;
        scaling[2] /= 100.0;

        FbxMesh* mesh = (FbxMesh*)attribute;
        const int vert_count = mesh->GetControlPointsCount();

        FbxVector4* vertices = mesh->GetControlPoints();
        FbxVector4* transformed_vertices = new FbxVector4[vert_count];

        FbxAMatrix transform(translation, rotation, scaling);

        // create polygons
        for (int i = 0; i < vert_count; i++)
        {
            transformed_vertices[i] = transform.MultT(vertices[i]);
            for (int j = 0; j < 3; ++j)
            {
                bounds[0][j] = FbxMin(bounds[0][j], transformed_vertices[i][j]);
                bounds[1][j] = FbxMax(bounds[1][j], transformed_vertices[i][j]);
            }
        }

        delete[] transformed_vertices;
    }
}

void rasterizeTriangles(FbxScene* scene, FbxVector4* bounds,
                        FbxVector4& sample_size, FbxVector4& inversed_sample_size,
                        ofstream& output_file)
{
#ifdef DEBUG_INFO
    cout << "Start processing fbx nodes ..." << endl;
#endif

    FbxNode* root_node = scene->GetRootNode();
    int child_count = root_node->GetChildCount();

#ifdef DEBUG_INFO
    cout << child_count << " nodes in total" << endl;
#endif

    unsigned int sample_count = 0;
    for (int index = 0; index < child_count; ++index)
    {
#ifdef DEBUG_INFO
        cout << "Start sampling " << index << "/" << child_count << " node" << endl;
#endif
        FbxNode* child_node = root_node->GetChild(index);
        if (child_node == nullptr) continue;
        FbxNodeAttribute* attribute = child_node->GetNodeAttribute();
        FbxNodeAttribute::EType attributeType = attribute->GetAttributeType();
        if (attributeType != FbxNodeAttribute::eMesh) continue;
        FbxDouble3 translation = child_node->LclTranslation.Get();
        FbxDouble3 rotation = child_node->LclRotation.Get();
        FbxDouble3 scaling = child_node->LclScaling.Get();

        translation[0] /= 100.0;
        translation[1] /= 100.0;
        translation[2] /= 100.0;

        scaling[0] /= 100.0;
        scaling[1] /= 100.0;
        scaling[2] /= 100.0;

        FbxMesh* mesh = (FbxMesh*)attribute;
        const int vert_count = mesh->GetControlPointsCount();

        FbxVector4* vertices = mesh->GetControlPoints();
        FbxVector4* transformed_vertices = new FbxVector4[vert_count];

        FbxAMatrix transform(translation, rotation, scaling);

        // create polygons
        for (int i = 0; i < vert_count; i++)
        {
            transformed_vertices[i] = transform.MultT(vertices[i]);
        }

        const int tri_count = mesh->GetPolygonCount();

#ifdef DEBUG_INFO
        cout << "(" << index << "/" << child_count << ") " << tri_count << " triangles in total" << endl;
#endif

        int* verts = mesh->GetPolygonVertices();
        for (int i = 0; i < tri_count; ++i)
        {
#ifdef DEBUG_INFO
            cout << "(" << index << "/" << child_count << ") " << "Sampling " << i << "/" << tri_count << " triangle: ";
#endif
            int poly_vert_start = mesh->GetPolygonVertexIndex(i);
            int* poly_verts = &verts[poly_vert_start];
            int vert_count = mesh->GetPolygonSize(i);
            assert(vert_count == 3);
            sample_count += rasterizeTriangle(transformed_vertices[poly_verts[0]],
                              transformed_vertices[poly_verts[1]],
                              transformed_vertices[poly_verts[2]],
                              bounds[0],
                              sample_size,
                              inversed_sample_size,
                              output_file);
        }

        delete[] transformed_vertices;
    }

#ifdef DEBUG_INFO
    cout << "Finished sampling: " << sample_count << " samples in total" << endl;
#endif
}

int main(int argc, char* argv[])
{
    // load fbx
    FbxManager* fbx_manager{nullptr};
    FbxScene* fbx_scene{nullptr};

    string output_file_name(argv[1]);
    ofstream output_file(output_file_name + ".bin", ios::binary);
    if (output_file.is_open() == false)
    {
        return -1;
    }

    initializeSdkObjects(fbx_manager, fbx_scene, "");

#ifdef DEBUG_INFO
    cout << "Start loading fbx file ..." << endl;
#endif

    loadScene(fbx_manager, fbx_scene, argv[1]);

#ifdef DEBUG_INFO
    cout << "Finished loading fbx file" << endl;

    cout << "Start calculating mesh bounding ..." << endl;
#endif

	FbxVector4 bounds[2];
	const double max_double = std::numeric_limits<double>::max();
	bounds[0] = FbxVector4(max_double, max_double, max_double);
	bounds[1] = FbxVector4(-max_double, -max_double, -max_double);
	calculateBounds(fbx_scene, bounds);

#ifdef DEBUG_INFO
    cout << "Finished calculating mesh bounding: ("
         << bounds[0][0] << ", " << bounds[0][1] << ", " << bounds[0][2] << ") ("
         << bounds[1][0] << ", " << bounds[1][1] << ", " << bounds[1][2] << ") " << endl;
#endif

    const FbxVector4 margin(0.5, 0.5, 0.5, 0.0);
    bounds[0] -= margin;
    bounds[1] += margin;
    FbxVector4 sample_size_vec(atof(argv[2]), atof(argv[3]), atof(argv[4]));
    FbxVector4 inversed_sample_size(1.0 / atof(argv[2]), 1.0 / atof(argv[3]), 1.0 / atof(argv[4]));

    float sample_size[3] = {0.f, 0.f, 0.f};
    int scene_bound[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
    {
        sample_size[i] = (float)sample_size_vec[i];
        scene_bound[i] = (int)::FbxFloor((bounds[1][i] - bounds[0][i]) * inversed_sample_size[i]) + 1;
    }

    output_file.write((char*)sample_size, sizeof(float) * 3);
    output_file.write((char*)scene_bound, sizeof(int) * 3);

    rasterizeTriangles(fbx_scene, bounds, sample_size_vec, inversed_sample_size, output_file);

    output_file.close();

    destroySdkObjects(fbx_manager, true);
}