#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{
    
    Vector2D lerp(Vector2D p0, Vector2D p1, float t)
    {
        Vector2D result;
        result.x = (1.0 - t) * p0.x + t * p1.x;
        result.y = (1.0 - t) * p0.y + t * p1.y;
        return result;
    }
    
    Vector3D Three_lerp(Vector3D p0, Vector3D p1, float t)
    {
        Vector3D result;
        result.x = (1.0 - t) * p0.x + t * p1.x;
        result.y = (1.0 - t) * p0.y + t * p1.y;
        result.z = (1.0 - t) * p0.z + t * p1.z;
        return result;
    }
    
    void BezierCurve::evaluateStep()
    {
        // TODO Part 1.
        // Perform one step of the Bezier curve's evaluation at t using de Casteljau's algorithm for subdivision.
        // Store all of the intermediate control points into the 2D vector evaluatedLevels.
        vector<Vector2D> points = evaluatedLevels[evaluatedLevels.size() - 1];
        vector<Vector2D> evaluated;
        if (points.size() <= 1) {
            return;
        }
        for (int i = 0; i < points.size() - 1; i++) {
            Vector2D p0 = points[i];
            Vector2D p1 = points[i + 1];
            evaluated.push_back(lerp(p0, p1, t));
        }
        evaluatedLevels.push_back(evaluated);
        return;
    }
    
    
    Vector3D BezierPatch::evaluate(double u, double v) const
    {
        // TODO Part 2.
        // Evaluate the Bezier surface at parameters (u, v) through 2D de Casteljau subdivision.
        // (i.e. Unlike Part 1 where we performed one subdivision level per call to evaluateStep, this function
        // should apply de Casteljau's algorithm until it computes the final, evaluated point on the surface)
        vector<Vector3D> q;
        for (int i = 0; i < controlPoints.size(); i++) {
            q.push_back(evaluate1D(controlPoints[i], u));
        }
        
        return evaluate1D(q, v);
    }
    
    Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> points, double t) const
    {
        // TODO Part 2.
        // Optional helper function that you might find useful to implement as an abstraction when implementing BezierPatch::evaluate.
        // Given an array of 4 points that lie on a single curve, evaluates the Bezier curve at parameter t using 1D de Casteljau subdivision.
        vector<Vector3D> Three_points;
        vector<Vector3D> Two_points;
        for (int i = 0; i < points.size() - 1; i++) {
            Three_points.push_back(Three_lerp(points[i], points[i + 1], t));
        }
        for (int i = 0; i < Three_points.size() - 1; i++) {
            Two_points.push_back(Three_lerp(Three_points[i], Three_points[i + 1], t));
        }
        return Three_lerp(Two_points[0], Two_points[1], t);
    }
    
    
    
    Vector3D Vertex::normal( void ) const
    {
        // TODO Part 3.
        // TODO Returns an approximate unit normal at this vertex, computed by
        // TODO taking the area-weighted average of the normals of neighboring
        // TODO triangles, then normalizing.
        Vector3D n(0,0,0);
        HalfedgeCIter h = halfedge();
        h = h->twin();
        HalfedgeCIter h_orig = h;
        do {
            Vector3D v0 = h->vertex()->position;
            Vector3D v1 = h->next()->twin()->vertex()->position;
            n += cross(v1, v0);
            h = h->next()->twin();
        } while (h != h_orig);
        return n.unit();
    }
    
    EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
    {
        // TODO Part 4.
        // TODO This method should flip the given edge and return an iterator to the flipped edge.
        HalfedgeIter h = e0->halfedge();
        if (h->isBoundary() || h->twin()->isBoundary()) {
            return e0;
        }
        
        if (h->vertex()->halfedge() == h) {
            h->vertex()->halfedge() = h->twin()->next();
        }
        if (h->twin()->vertex()->halfedge() == h->twin()) {
            h->twin()->vertex()->halfedge() = h->next();
        }
        
        //Change source Vertex
        h->vertex() = h->next()->next()->vertex();
        h->twin()->vertex() = h->twin()->next()->next()->vertex();
        
        //Reassign next for all half edges
        HalfedgeIter t0 = h->next()->next();
        HalfedgeIter t1 = h->twin()->next()->next();
        h->next()->next()->next() = h->twin()->next();
        h->twin()->next()->next()->next() = h->next();
        h->next()->next() = h;
        h->twin()->next()->next() = h->twin();
        h->twin()->next() = t0;
        h->next() = t1;
        h->twin()->twin() = h;
        
        //Reassign faces/faces' halfedges
        h->face()->halfedge() = h;
        h->next()->face() = h->face();
        h->next()->next()->face() = h->face();
        h->twin()->face()->halfedge() = h->twin();
        h->twin()->next()->face()= h->twin()->face();
        h->twin()->next()->next()->face()= h->twin()->face();
        
        h->edge()->halfedge() = h;
        h->next()->edge()->halfedge() = h->next();
        h->next()->next()->edge()->halfedge() = h->next()->next();
        h->twin()->next()->edge()->halfedge() = h->twin()->next();
        h->twin()->next()->next()->edge()->halfedge() = h->twin()->next()->next();
        
        h->next()->next()->setNeighbors(h, h->next()->next()->twin(), h->next()->next()->vertex(), h->next()->next()->edge(), h->next()->next()->face());
        h->next()->setNeighbors(h->next()->next(), h->next()->twin(), h->next()->vertex(), h->next()->edge(), h->next()->face());
        h->setNeighbors(h->next(), h->twin(), h->vertex(), h->edge(), h->face());
        h->twin()->next()->next()->setNeighbors(h->twin(), h->twin()->next()->next()->twin(), h->twin()->next()->next()->vertex(), h->twin()->next()->next()->edge(), h->twin()->next()->next()->face());
        h->twin()->next()->setNeighbors(h->twin()->next()->next(), h->twin()->next()->twin(), h->twin()->next()->vertex(), h->twin()->next()->edge(), h->twin()->next()->face());
        h->twin()->setNeighbors(h->twin()->next(), h, h->twin()->vertex(), h->twin()->edge(), h->twin()->face());
        
        return e0;
    }
    
    VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
    {
        // TODO Part 5.
        // TODO This method should split the given edge and return an iterator to the newly inserted vertex.
        // TODO The halfedge of this vertex should point along the edge that was split, rather than the new edges.
        HalfedgeIter h = e0->halfedge();
        if (h->isBoundary() || h->twin()->isBoundary()) {
            return e0->halfedge()->vertex();
        }
        //Create new vertex
        VertexIter v = newVertex();
        Vector3D midpoint = (h->vertex()->position + h->twin()->vertex()->position)/2.0;
        v->position = midpoint;
        
        //Create new halfedges and reassign twins
        //dcm
        HalfedgeIter md = newHalfedge();
        HalfedgeIter dm = newHalfedge();
        HalfedgeIter mb = newHalfedge();
        HalfedgeIter bm = newHalfedge();
        HalfedgeIter ma = newHalfedge();
        HalfedgeIter am = newHalfedge();
        md->twin() = dm;
        dm->twin() = md;
        mb->twin() = bm;
        bm->twin() = mb;
        ma->twin() = am;
        am->twin() = ma;
        
        //Assign vertices
        dm->vertex() = h->next()->next()->vertex();
        md->vertex() = v;
        bm->vertex() = h->twin()->vertex();
        mb->vertex() = v;
        am->vertex() = h->twin()->next()->next()->vertex();
        ma->vertex() = v;
        h->twin()->vertex()->halfedge() = bm;
        h->vertex()->halfedge() = h->twin()->next();
        h->twin()->vertex() = v;
        
        bm->vertex()->halfedge() = bm;
        h->vertex()->halfedge() = h;
        
        //Reassign next
        md->next() = h->next()->next();
        mb->next() = h->next();
        ma->next() = h->twin()->next()->next();
        
        mb->next()->next() = dm;
        dm->next() = mb;
        h->next() = md;
        
        bm->next() = ma;
        h->twin()->next()->next()->next() = bm;
        am->next() = h->twin();
        h->twin()->next()->next() = am;
        
        //Reassign Edges
        EdgeIter eAM = newEdge();
        EdgeIter eBM = newEdge();
        EdgeIter eDM = newEdge();
        
        eAM->halfedge() = am;
        eBM->halfedge() = bm;
        eDM->halfedge() = dm;
        
        am->edge() = eAM;
        ma->edge() = eAM;
        bm->edge() = eBM;
        mb->edge() = eBM;
        md->edge() = eDM;
        dm->edge() = eDM;
        h->edge()->halfedge() = h;
        h->twin()->edge() = h->edge();
        
        //Reassign Faces
        FaceIter mbd = newFace();
        FaceIter mab = newFace();
        
        h->face()->halfedge() = md;
        h->next()->face() = h->face();
        h->next()->next()->face() = h->face();
        h->twin()->face()->halfedge() = h->twin();
        h->twin()->next()->face() = h->twin()->face();
        h->twin()->next()->next()->face() = h->twin()->face();
        
        mbd->halfedge() = mb;
        mab->halfedge() = ma;
        mb->face() = mbd;
        mb->next()->face() = mbd;
        mb->next()->next()->face() = mbd;
        ma->face() = mab;
        ma->next()->face() = mab;
        ma->next()->next()->face() = mab;
        
        //Set halfEdge for Vertex
        v->isNew = true;
        eAM->isNew = true;
        eDM->isNew = true;
        eBM->isNew = false;
        h->edge()->isNew = false;
        
        v->halfedge() = h->twin();
        
        return v;
    }
    
    
    
    void MeshResampler::upsample( HalfedgeMesh& mesh )
    {
        // TODO Part 6.
        // This routine should increase the number of triangles in the mesh using Loop subdivision.
        // Each vertex and edge of the original surface can be associated with a vertex in the new (subdivided) surface.
        // Therefore, our strategy for computing the subdivided vertex locations is to *first* compute the new positions
        // using the connectity of the original (coarse) mesh; navigating this mesh will be much easier than navigating
        // the new subdivided (fine) mesh, which has more elements to traverse. We will then assign vertex positions in
        // the new mesh based on the values we computed for the original mesh.
        
        // TODO Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
        // TODO and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
        // TODO a vertex of the original mesh.
        for(VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
            v->isNew = false;
            
            HalfedgeCIter h = v->halfedge();
            float neighbors = 0.0;
            Vector3D neighbor_position_sum(0.0, 0.0, 0.0);
            Vector3D orig_position = v->position;
            do {
                HalfedgeCIter h_twin = h->twin();
                VertexCIter v0 = h_twin->vertex();
                neighbor_position_sum += v0->position;
                h = h_twin->next();
                neighbors = 1.0 + neighbors;
            } while(h != v->halfedge());
            
            float u = 3.0/(8.0 * (float) neighbors);
            
            if (neighbors == 3.0) {
                u = 3.0/16.0;
            }
            
            Vector3D newPos = (1.0 - ((float) neighbors * u)) * orig_position + u * neighbor_position_sum;
            v->newPosition = newPos;
        }
        
        // TODO Next, compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
        for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
            e->isNew = false;
            HalfedgeIter h = e->halfedge();
            Vector3D ab = h->vertex()->position + h->twin()->vertex()->position;
            Vector3D cd = h->next()->next()->vertex()->position + h->twin()->next()->next()->vertex()->position;
            
            Vector3D newPos = (3.0/8.0) * ab + (1.0/8.0) * cd;
            e->newPosition = newPos;
        }
        
        // TODO Next, we're going to split every edge in the mesh, in any order.  For future
        // TODO reference, we're also going to store some information about which subdivided
        // TODO edges come from splitting an edge in the original mesh, and which edges are new,
        // TODO by setting the flat Edge::isNew.  Note that in this loop, we only want to iterate
        // TODO over edges of the original mesh---otherwise, we'll end up splitting edges that we
        // TODO just split (and the loop will never end!)
        
        int numEdges = 0;
        for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
            numEdges += 1;
        }
        EdgeIter e = mesh.edgesBegin();
        for (int i = 0; i < numEdges + 1; i++) {
            EdgeIter nextEdge = e;
            nextEdge++;
            
            if (!e->isNew) {
                mesh.splitEdge(e);
            }
            e = nextEdge;
        }
        
        // TODO Now flip any new edge that connects an old and new vertex.
        e = mesh.edgesBegin();
        while (e != mesh.edgesEnd()) {
            EdgeIter nextEdge = e;
            nextEdge++;
            
            HalfedgeCIter h = e->halfedge();
            bool newVertex = (h->vertex()->isNew && !h->twin()->vertex()->isNew) || (!h->vertex()->isNew && h->twin()->vertex()->isNew);
            if (h->edge()->isNew && newVertex) {
                mesh.flipEdge(e);
            }
            e = nextEdge;
        }
        
        // TODO Finally, copy the new vertex positions into final Vertex::position.
        for(VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
            if (!v->isNew) {
                v->position = v->newPosition;
            } else {
                HalfedgeCIter h = v->halfedge();
                v->position = h->edge()->newPosition;
            }
        }
        
        return;
    }
}
