def normal(vec1,vec2):
    norv = [vec1[1]*vec2[2]-vec1[2]*vec2[1],vec1[2]*vec2[0]-vec1[0]*vec2[2],vec1[0]*vec2[1]-vec1[1]*vec2[0]]
    return norv

def angle(vec1,vec2):
    import math
    from numpy import radians,sin,cos,array,dot,sqrt
    cosa=dot(vec1,vec2)/sqrt(dot(vec1,vec1)*dot(vec2,vec2))
    return math.degrees(math.acos(cosa))

def rotation_zxz(phi,theta,psi):
    
    #phi: angle between lab frame x and molecular frame a
    #theta: angle between lab frame z and molecular frame c
    #psi: angle between the crossing line of plane xy and plane ab and molecular frame a
    
    from numpy import radians,sin,cos,array
    psi = radians(psi)
    phi = radians(phi)
    theta = radians(theta)

    a11 = cos(psi)*cos(phi) - cos(theta)*sin(psi)*sin(phi)
    a12 = -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi)
    a13 = sin(theta)*sin(phi)
    a21 = cos(psi)*sin(phi) + cos(theta)*sin(psi)*cos(phi)
    a22 = -sin(psi)*sin(phi) + cos(theta)*cos(psi)*cos(phi)
    a23 = -sin(theta)*cos(phi)
    a31 = sin(theta)*sin(psi)
    a32 = sin(theta)*cos(psi)
    a33 = cos(theta)

    rotmat = array([[ a11, a12, a13 ], [a21, a22, a23 ],[ a31, a32, a33 ]], dtype=float)
    return rotmat

def rotation_zyz(phi,theta,psi):
    #phi: angle between lab frame y and molecular frame b
    #theta: angle between lab frame z and molecular frame c
    #psi: angle between the crossing line of plane xy and plane ab and molecular frame a

    from numpy import radians,sin,cos,array
    phi = radians(phi)
    theta = radians(theta)
    psi = radians(psi)

    a11 = - sin(psi)*sin(phi) + cos(theta)*cos(psi)*cos(phi)
    a12 = -cos(psi)*sin(phi) - cos(theta)*sin(psi)*cos(phi)
    a13 = sin(theta)*cos(phi)
    a21 = sin(psi)*cos(phi) + cos(theta)*cos(psi)*sin(phi)
    a22 = cos(psi)*cos(phi) - cos(theta)*sin(psi)*sin(phi)
    a23 = sin(theta)*sin(phi)
    a31 = -sin(theta)*cos(psi)
    a32 = sin(theta)*sin(psi)
    a33 = cos(theta)

    rotmat = array([[ a11, a12, a13 ],
                    [ a21, a22, a23 ],
                    [ a31, a32, a33 ]], dtype=float)
    return rotmat

def plane(Obj,p0,p1,p2):

    from numpy import array

    Obj.translate_coordinates(-Obj.center_of_mass)
    vec1 = p1-p0 # x-axis
    vec2 = p2-p0 # another vector in the plane
    nor = normal(vec1,vec2)
    px = array([1,0,0])
    pz = array([0,0,1])
    share = normal(nor,pz)

    gamma = angle(vec1,share) 
    beta = angle(nor,pz)
    alpha = angle(share,px)

    rotmax = rotation_zxz(alpha,beta,gamma)
    Obj.rotate_coordinates(rotmax)
    Obj.writeCoords()
