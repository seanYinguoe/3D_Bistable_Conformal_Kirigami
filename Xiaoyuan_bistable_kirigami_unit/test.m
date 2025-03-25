
import Rhino.Geometry as rg
import re
from math import sqrt, pi, tan,cos,sin,degrees, atan2, atan
import rhinoscriptsyntax as rs
from collections import Counter

import os
import csv

import ghpythonlib.treehelpers as th
import copy

REX_VT = re.compile(r"vt\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
REX_FACE = re.compile(r"f\s+(\d+)\/(\d+)\s+.*?(\d+)\/(\d+)\s+.*?(\d+)\/(\d+).*?")

def vecCross(a, b):
    if len(a)==3:
        c = [a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0]]
    elif len(a)==2:
        c = [0, 0, a[0]*b[1] - a[1]*b[0]]
    return c

def norm(a):
    return sqrt(sum([a[i]**2 for i in [0,1,2]]))
     
def vecSub(a,b): # a - b
    if len(a)==3:
        return [a[i]-b[i] for i in [0,1,2]]
    elif len(a)==2:
        return [a[i]-b[i] for i in [0,1]]
            
def triangle_areas(v, f):
    areas = []
    for i,j,k in f:
        area2 = norm(vecCross(vecSub(v[i], v[j]), vecSub(v[i], v[k])))
        areas.append(area2/2.)
    return areas

def get_floppy_triangles(f_grid):
    '''
    filter the triangles that has at least one vertex that is not connected to any other triangles
    '''
    f_grid_flat = [item for sublist in f_grid for item in sublist]
    
    res= Counter (f_grid_flat)
    
    idE = []
    
    for key, item in res.items():  # for name, age in dictionary.iteritems():  (for Python 2.x)
        if item == 1:
            idE.append(key)
            
    return idE

def removeId(f_grid, idE):
    '''
    Remove the floppy triangles by replacing them with [-1] so that the numbering doesn't change
    '''
    temp = []
    id_y = []
    id_n = []
    for i,fi in enumerate(f_grid):
        if fi != [-1]: # if it is not already filtered out
            if fi[0] not in idE and fi[1] not in idE and fi[2] not in idE: #check if any of each face has a floppy vertex
                temp.append(fi)
                id_y.append(i)
            else: # if there is one, then change the face to [-1]
                temp.append([-1])
                id_n.append(i)
        else:
            temp.append([-1])
    return id_y,id_n, temp

def normList(l):
    return [li/sum(l) for li in l]
    
def distClosest(needle, hay, d):
    '''
    return all the indices in the list that is close to a value by distance d (1-D)
    '''
    bL = [abs(x-needle)<d for x in hay]
    return [i for i, x in enumerate(bL) if x]
    
def ppClosest(h, n):
    '''
    Find the point in h that n is closest to
    '''
    tempD = [rs.Distance(hi,n) for hi in h]
    if min(tempD)>1:
        errorTog = True
    else:
        errorTog = False
    return tempD.index(min(tempD)), errorTog

def noneList(x,a):
    ''' change the indices a of x to [-1] '''
    for ai in a:
        x[ai] = [-1]
        
def noneListR(a):
    ''' remove [-1] '''
    idxD = dict()
    idx = 0
    temp = []
    for i, ai in enumerate(a):
        if ai != [-1]:
            temp.append(ai)
            idxD[i]=idx
            idx+=1
            
    return temp, idxD
        
    
def argmax(l):
    return l.index(max(l))
    
def argmin(l):
    return l.index(min(l))
    
def remove_unused_verts(v, f):
    """
    Takes a mesh in form of verts and triangles and
    returns a new mesh with the verts removed that are not
    referenced in the triangles.
    """
    
    v_new = []
    f_new = []
    c_new = []
    
    for fi in f:
        (i0, i1, i2) = fi
        
        # Round is cheap method to check do "is approximately in"
        v1 = [round(x,8) for x in v[i0]]
        v2 = [round(x,8) for x in v[i1]]
        v3 = [round(x,8) for x in v[i2]]
        
        # if v1,2,3 is not in v_new then it is a new point that should be appended
        if v1 not in v_new:
            v_new.append(v1)
        if v2 not in v_new:
            v_new.append(v2)
        if v3 not in v_new:
            v_new.append(v3)
        
        # faces then append the id of all three points, 
        # now that they have to be in v_new
        f_new.append([v_new.index(v1), v_new.index(v2), v_new.index(v3)])
        
    return v_new, f_new

def load_obj(obj_path):
    file = open(obj_path)
    lines = file.readlines()

    f = []
    v = []
    uv = []
    
    t_mesh  = rg.Mesh()  # Target mesh
    uv_mesh = rg.Mesh()  # UV (flat) mesh
    
    f_lines  = [line for line in lines if line.startswith("f")]
    v_lines  = [line for line in lines if (line.startswith("v") and line.find("n")==-1 and line.find("t")==-1)]
    uv_lines = [line for line in lines if line.startswith("vt")]
    
    for line in v_lines:
        v.append([float((line.split(' '))[1]),float((line.split(' '))[2]),float((line.split(' '))[3])])
        t_mesh.Vertices.Add(rg.Point3d(*v[-1]))
    
    uvs = [[float(coord) for coord in REX_VT.findall(line)[0]] for line in uv_lines]
    uv_indices = [None] * len(uvs)
    
    for line in f_lines:
            f.append([int(line.split(' ')[1].split('/')[0])-1,int(line.split(' ')[2].split('/')[0])-1,int(line.split(' ')[3].split('/')[0])-1])
            
            mf = rg.MeshFace(*f[-1])
            
            uv_mesh.Faces.AddFace(mf) # the two meshes have the same connectivities
            t_mesh.Faces.AddFace(mf)  # the two meshes have the same connectivities
            
            res_str = REX_FACE.findall(line)[0]
            res_int = [int(index) - 1 for index in res_str]

            v1, uv1, v2, uv2, v3, uv3 = res_int

            try:
                uv_indices[v1] = uv1
                uv_indices[v2] = uv2
                uv_indices[v3] = uv3
            except Exception as e:
                print("uv1", uv1, "uv2", uv2, "uv3", uv3)
                print(len(uvs))
                raise e

    uvs_sorted = [uvs[i] for i in uv_indices]
    
    [uv_mesh.Vertices.Add(rg.Point3d(uvs[0],uvs[1],0)) for uvs in uvs_sorted]
    uv_mesh.Compact()
    t_mesh.Normals.ComputeNormals()
    t_mesh.Compact()
    
    return [v,f,uvs_sorted,uv_mesh,t_mesh]

def load_mesh_and_scale_facs(parametrized_obj_path, verbose):
    [v,f,uvs,uv_mesh,t_mesh]=load_obj(parametrized_obj_path)

    tri_areas3d = triangle_areas(v, f)
    tri_areas2d = triangle_areas(uvs, f)

    scale_facs = [x/y for x,y in zip(tri_areas3d, tri_areas2d)]

    zero_area_2dtris = [tri_areas2d < 1e-8 for x in tri_areas2d]

    for i,x in enumerate(zero_area_2dtris):
        if x == True:
            scale_facs[i] = None

    if sum(zero_area_2dtris) > 0:
        warnings.warn(
            "Something is wrong with your UV-map. Some triangles have near-zero or zero area."
            "This leads to infinite scale factors."
        )

    if verbose:
        print("Scale facs min: "+str(min(scale_facs)) +", max: "+str(max(scale_facs)))

    return v, f, uvs, scale_facs, uv_mesh, t_mesh

def flip(i):
    if i==0:
        return 1
    elif i==1:
        return 0

# in even last one, i dont exclude the upper triangle, which causes issues
def generate_regular_tri_grid(rows, cols, l=1, ox=0, oy=0):
    v = []
    f = []
    c = []
    i_grid = []
    x_grid = []
    
    h = sqrt(3)/2*l # Height of the triangle
    
    idE = 0 # Counter for elements
    for row in range(rows):
        for col in range(cols):
            if row%2==0:
                v.append([l*col+ox,h*row+oy, 0])
            else:
                v.append([l*col+l/2+ox,h*row+oy, 0])
                
    for row in range(rows-1):
        iList = []
        for col in range(cols-1):
            
            if row%2==0:
                m = 1
                n = 0
                iList.append(0)
                iList.append(1)
                
            else:
                m = 0
                n = 1
                iList.append(1)
                iList.append(0)
        
            n0 = col + row*cols
            n1 = col + row*cols + 1
            n2 = col + (row+1)*cols + n
            # iso triangle #1
            fL1 = [n0,n1,n2]
        
            c0 = (v[n0][0]+v[n1][0])/2
            c1 = v[n0][1] + sqrt(3)/6*l
            # center of tri #1
            cL1 = [c0,c1,0]
            
            n0 = col + row*cols + m
            n1 = col + (row+1)*cols + 1
            n2 = col + (row+1)*cols
            
            # iso triangle #2
            fL2 = [n0,n1,n2]
            
            c0 = (v[n1][0]+v[n2][0])/2
            c1 = v[n1][1] - sqrt(3)/6*l
            
            # center of tri #2
            cL2 = [c0,c1,0]

            # this generation always starts with an a triangle with bottom edge aligned with the xaxis first
            if row%2==0:
                    
                f.append(fL1)
                c.append(cL1)
                f.append(fL2)
                c.append(cL2)
                
                # Find the ID of the neighbouring elements
                if row == 0:
                    if col == 0:
                        x_grid.append([idE+1])
                    else:
                        x_grid.append([idE-1, idE+1])
                else:
                    if col ==0:
                        x_grid.append([idE+1, idE-(cols-1)*2])
                    else:
                        x_grid.append([idE-1, idE-(cols-1)*2, idE+1])
                     
                idE+=1   
                
                # Find the ID of the neighbouring elements
                if row == rows-2:
                    if col == cols-2:
                        x_grid.append([idE-1])
                    else:
                        x_grid.append([idE-1, idE+1])
                else:
                    if col == cols-2:
                        x_grid.append([idE-1, idE+(cols-1)*2])
                    else:
                        x_grid.append([idE-1, idE+(cols-1)*2, idE+1])
                        
                idE+=1
                
            else:
                f.append(fL2)
                c.append(cL2)
                
                f.append(fL1)
                c.append(cL1)
                
                # Find the ID of the neighbouring elements
                if row == rows - 2:
                    if col ==0:
                        x_grid.append([idE+1])
                    else:
                        x_grid.append([idE-1, idE+1])

                else:
                    if col ==0:
                        x_grid.append([idE+1, idE+(cols-1)*2])
                    else:
                        x_grid.append([idE-1, idE+(cols-1)*2, idE+1])
                    
                idE+=1
                
                # Find the ID of the neighbouring elements
                if row == 0:
                    # But this wont happen unless we flip the order or even and odd
                    if col == cols - 2:
                        x_grid.append([idE-1])
                    else:
                        x_grid.append([idE-1, idE+1])
                else:
                    if col == cols-2:
                        x_grid.append([idE-1, idE-(cols-1)*2])
                    else:
                        x_grid.append([idE-1, idE-(cols-1)*2, idE+1])
                        
                idE+=1
                
        i_grid.append(iList)
        
    i_grid = [item for sublist in i_grid for item in sublist]
    
    return v,f,c,i_grid, x_grid

def generate_overlay_grid(rsv2d, params):
    """
    Takes a flattened mesh and returns the verts, faces
    for a triangular grid that envelops the input mesh.
    """
    boundv2d = rs.PointArrayBoundingBox(rsv2d)

    cols = int(round((boundv2d[2][0] - boundv2d[0][0])/in_l_o,0)+4)
    rows = int(round((boundv2d[2][1] - boundv2d[0][1])/(sqrt(3)/2*in_l_o),0)+4)

    ox  = boundv2d[0][0] - in_l_o * 2 + params['grid_offset'][0]
    oy  = boundv2d[0][1] - sqrt(3)/2 * in_l_o * 2 + params['grid_offset'][1]
 
    if rows + cols < 5:
        warnings.warn("The design is very small. Check your model size and units (mm assumed).")

    v_grid, f_grid, c_grid, i_grid, x_grid = generate_regular_tri_grid(rows, cols, in_l_o, ox, oy)
    
    return v_grid, f_grid, c_grid, i_grid, x_grid

def readData(fName):
    # All this stuff is only done once, at the time of module loading

    with open (fName) as csvfile:
        temp = csv.reader(csvfile, delimiter=',')
        rows = []
        for row in temp:
            rows.append(row)
            
    ts          = [float(x) for x in rows[0]]
    thetas      = [float(x) for x in rows[1]]
    expansions  = [float(x) for x in rows[2]]
    energies    = [float(x) for x in rows[3]]
    stiffnesses = [float(x) for x in rows[4]]
    
    # output the min and max of the parameter space
    print('ts_norm (t/l) min = ' + str(min(ts))     + ', max = ' + str(max(ts)))
    print('thetas  (rad) min = ' + str(min(thetas)) + ', max = ' + str(max(thetas)))
   
    return [ts, thetas, expansions, energies, stiffnesses]

def get_most_stable_for_expansion(target_expansion, data, tol_expansion=0.005, mode='stiffness', verbose=False, debug=False):

    ts          = data[0]
    thetas      = data[1]
    expansions  = data[2]
    energies    = data[3]
    stiffnesses = data[4]
    
    # Get the id of neighbours for each target exp within the tol
    exp_neigh   = [distClosest(x, expansions, tol_expansion) for x in target_expansion]
    
    id_selected = []
    select_neigh = []
    for ids in exp_neigh:
        # If none then not bistable
        if not ids:
            select_neigh.append(None)
            return
        if mode =='energy':
            # get barrier height of these close neighbours
            barrier_neigh = [energies[idi] for idi in ids]
            select_neigh  = barrier_neigh
            
        elif mode =='stiffness':
            # get barrier height of these close neighbours
            stiff_neigh   = [stiffnesses[idi] for idi in ids]
            select_neigh  = stiff_neigh
        
        id_selected.append(ids[argmax(select_neigh)])
        
    return id_selected
    
def interp_and_rescale(v_grid, f_grid, c_grid, x_grid, uv_mesh, i_grid, scale_facs, scale_min, scale_max, geomType, dbName, mode='energy',verbose=True):
    """
    Takes a regular triangular grid and a flattened triangular mesh
    and interpolates the scale factors from the mesh to the regular grid.

    for_expanded_cells: Shift the scale factors such that we get parameters for opened cells that close
                        instead of closed ones that open 

    use_bb: Use the bounding box of the uv map instead of the uv map itself to fit the tri grid into.
            Useful for cylinders, i.e. if you need the top and bottom of the pattern to be gluable together. Will
            introduce some error however, it's better to adapt the boundary in BFF if possible.
    """

    # centroids flat is the coordinates of the center of the flattened triangles
    # scale_facs is calculated per uv triangle [num mesh elem x 1]
    centroids_uv = rs.MeshFaceCenters(uv_mesh)

    # if not using the bounding box, it will use the uv map instead.
    # the FieldSampler checks for holes or concave parts of the uv map
    # if the regular triangles are inside the mesh then True, otherwise false.
    
    inside_uv_tris = [uv_mesh.ClosestMeshPoint(rs.coerce3dpoint(rs.AddPoint(x)),0.0000000001) for x in c_grid]
    
    inside_tris = [i for i,tri in enumerate(inside_uv_tris) if tri is None]
    
    # replace all the triangles not inside by -1
    noneList(f_grid,inside_tris)
    noneList(c_grid,inside_tris)
    noneList(i_grid,inside_tris)
    noneList(x_grid,inside_tris)
    
    # get triangles with at least one vertex that's not shared by any other triangle
    idE = get_floppy_triangles(f_grid)
    
    # remove those triangles from the connectivity by -1
    idx, idxn, f_grid = removeId(f_grid, idE)
    
    # replace all the floppy tri by -1
    noneList(c_grid, idxn)
    noneList(i_grid, idxn)
    noneList(x_grid, idxn)
    
    # Actually remove all the -1, and keep the coressponding indices
    c_grid,_    = noneListR(c_grid)
    i_grid,_    = noneListR(i_grid)
    x_grid,_    = noneListR(x_grid)
    f_grid,idxD = noneListR(f_grid)
    
    # do the same for x_grid, since it is a neighbouring tri relationship
    
    temp = []
    for xi in x_grid:
        tempL = [idxD[xii] for xii in xi if xii in idxD]
        temp.append(tempL)
    x_grid = temp
    # now calculate the scale factor by averaging the scale factors of the neighbouring mesh elements in the uv mesh
    
    scale_area = []
    
    for ci in c_grid:
        index = rs.PointCloudKNeighbors(centroids_uv, ci, 5)
        d = normList([1./rs.Distance(ci, centroids_uv[i]) for i in index[0]])
        scale_area.append(sum([scale_facs[idi]*di for idi,di in zip(index[0],d)]))
        
    # get length from area change by sqrting
    scale_length   = [sqrt(x) for x in scale_area]

    if verbose == True:
        print("Min factor (uv): " + str(min(scale_length)))
        print("Max factor (uv): " + str(max(scale_length)))
        
    # now shift the scale range such that the the largest scale is the maximum 
    # (this is a heuristic, it means that the flat shape will be the smallest)
    uv_scale = scale_max / max(scale_length)
    
    scale_length = [x*uv_scale for x in scale_length]

    if verbose ==True:
        print("Scaled by " + str(uv_scale) + " to")
        print("Scaled min factor (uv): " + str(min(scale_length)))
        print("Scaled max factor (uv): " + str(max(scale_length)))
        print ('')
        
    errorTog = False
    # this only happens with expanded cells
    if max(scale_length) > scale_max*1.05:
        print("Design exceeds maximal length scaling")
        errorTog = True
    if min(scale_length) < scale_min/1.05:
        print("Design exceeds minimum length scaling")
        errorTog = True

    if geomType == 'bass':
        param_out = []     # params
        
        data = readData(dbName)
        
        id_selected = get_most_stable_for_expansion(scale_length, data, mode=mode)
        if id_selected != None:
            for idi in id_selected:
                param_out.append([data[i][idi] for i in range(0,5)])
            
    elif geomType == 'aux':
        scale_area = [(s/max(scale_length))*2 for s in scale_length]
        param_out    = [auxPhi(x) for x in scale_area]
        
    v_out     = v_grid # coordinates
    f_out     = f_grid # connectivity
    c_out     = c_grid
    i_out     = i_grid
    x_out     = x_grid
    
    v_out, f_out = remove_unused_verts(v_out, f_out)
    
    return v_out, f_out, c_out, i_out, x_out, param_out, uv_scale, scale_length, centroids_uv, errorTog
   
def auxPhi (xf):
    if xf >= 2:
        phi = 0
    elif xf <= 1:
        phi = pi/3
    else:
        phi = atan2(sqrt(3*xf - 3*sqrt(-xf**4 + 4*xf**3 - 3*xf**2))/2 , -sqrt(3*xf - 3*sqrt(-xf**4 + 4*xf**3 - 3*xf**2))*(xf**2 + 2*sqrt(-xf**4 + 4*xf**3 - 3*xf**2))*sqrt(3)/(6*xf*(xf - 2)))
    return phi
 
def scaleRange(geomType):
    if geomType =='bass':
        return [1.15,1.75]
    if geomType =='aux':
        return [1.,2.]

scaleRange = scaleRange(geomType)

if not grid_offset:
    grid_offset = [0,0,0]
    
params = {
    "modelPath":     modelPath,
    "modelName":     modelName,
    "DBpath":        DBpath,
    "scale_minimal": scaleRange[0],
    "scale_maximal": scaleRange[1],
    "grid_offset":   [grid_offset[i] for i in [0,1]],      # Manual control for centering
    "selectMode":    selectMode,    #stiffness, energy
    "geomType":      geomType
}

# The UV map (called 2D) is scaled here prior to calculating the necessary global scaling
v, f, uv, scale_facs, uv_mesh, t_mesh = load_mesh_and_scale_facs(params['modelPath'] + '/'+params['modelName']+'/'+params['modelName']+"_flat.obj", verbose=False) #0.9725

rs_v  = rs.AddPoints(v)
rs_uv = rs.AddPoints(uv)
# Generate a regular triangular grid that is centered around the flattened uv mesh
v_grid, f_grid, c_grid, i_grid, x_grid = generate_overlay_grid(rs_uv, params)

rsv_grid = rs.AddPoints(v_grid)
rsc_grid = rs.AddPoints(c_grid)

[v_out, f_out, c_out, i_out, x_out, param_out, uv_scale, scale_length,centroids_uv,errorTog] = interp_and_rescale(
        v_grid,f_grid,c_grid,x_grid, uv_mesh,
        i_grid     = i_grid,
        scale_facs = scale_facs,
        scale_min  = params["scale_minimal"],
        scale_max  = params["scale_maximal"],
        geomType   = params["geomType"],
        dbName     = params['DBpath'],
        mode       = params["selectMode"])

rsv_grid = rs.AddPoints(v_out)
rsf_grid = [rs.AddPolyline([rsv_grid[f[0]],rsv_grid[f[1]],rsv_grid[f[2]],rsv_grid[f[0]]]) for f in f_out]
rsc_grid = rs.AddPoints(c_out)

paramsModel = params.items()

uv_mesh_edge=uv_mesh.GetNakedEdges()
c_out = rs.AddPoints(c_out)
[paramsModel,param_out,v_out,f_out,x_out] = [th.list_to_tree(x) for x in [paramsModel,param_out,v_out,f_out,x_out]]

if errorTog == True:
    msg = 'Design exceeds length scaling for bistable auxetics, try using auxetics instead.'
else:
    msg = 'No errors'