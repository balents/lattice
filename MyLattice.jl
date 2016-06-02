# Lattice generation

module MyLattice

export Lattice, coord, flavor, nns, nnns, num_sites, num_uc, adj
export bravais_full


# This will not be accessible outside the module
function bravais_int(Xvecs)
    # Xvecs should be a d by d matrix of integers giving winding vectors as columns
    # This function returns a unique set of integer coordinates of points defining a cluster given the integer winding numbers  

    const nmax = 40
    dim = size(Xvecs,1)
    if dim == 3
    	newxvecs = Xvecs
    elseif dim == 2
    	newxvecs = zeros(Int64,3,3)
    	newxvecs[1:2,1:2] = Xvecs
    	newxvecs[3,3] = 1
    else
        error("Lattice must have dimension 2 or 3")
    end
    tp = dot(newxvecs[:,1],cross(newxvecs[:,2],newxvecs[:,3]))

    numpts = abs(tp)
    function reduced_pt(pt)
        # gives a unique integer point inside a supercell defined by Xvecs
        # Reciprocal lattice for the supercell
        
        dim = length(pt)
        if dim == 2
            pt3d = [pt; 0]
        else
            pt3d = pt
        end 
        recip = Array{Float64}(3,3)
        for a in 1:3
            recip[:,a] = cross(newxvecs[:,mod(a,3)+1],newxvecs[:,mod(a+1,3)+1])/tp
        end
        fracpt = mod(round(recip'*float(pt3d),3),1)
        newpt = round(Int64,newxvecs*fracpt)[1:dim]
        return newpt
    end

    myset=Set{Vector{Int64}}()
    skip = false
    if dim == 3
        for i in 0:nmax
            for j in 0:nmax
                for k in 0:nmax
                    thispt = [i,j,k]              
                    modpt = reduced_pt(thispt)
                    push!(myset,modpt)
                    length(myset)==numpts && (skip = true; break)
                end
                skip && break
            end
            skip && break
        end
    else
    	for i in 0:nmax
    		for j in 0:nmax
    			thispt = [i,j]
    			modpt = reduced_pt(thispt)
    			push!(myset,modpt)
    			length(myset)==numpts && (skip = true; break)
    		end
    		skip && break
    	end
    end

    intpts = Array{Int64}(dim,numpts)
    for (i,vec) in enumerate(myset)
        intpts[:,i] = vec
    end

    return intpts
end

# Lattice type.  Stores all sorts of good stuff

type Lattice
    sites::Matrix{Float64}
    flavor::Vector{Int64}
    firstneighbors::Matrix{Int64}
    secondneighbors::Matrix{Int64}
    adjacency::Matrix{Int64}
    nsites::Int64
    nuc::Int64
    numnns::Int64
    numnnns::Int64
end

coord(lat::Lattice,i) = lat.sites[:,i]
flavor(lat::Lattice,i) = lat.flavor[i]
nns(lat::Lattice,i) = lat.firstneighbors[:,i]
nnns(lat::Lattice,i) = lat.secondneighbors[:,i]
num_sites(lat::Lattice) = lat.nsites
num_uc(lat::Lattice) = lat.nuc
adj(lat::Lattice,i,j) = lat.adjacency[i,j]

function bravais_full(avecs,basis,Xvecs)
    # Will construct a Lattice object describing a cluster of lattice defined by these parameters
    # avecs should be a d by d matrix of primitive lattice vectors arranged in columns
    # basis should be a d by n matrix of basis vectors for n sites in a unit cell.
    # Xvecs should be a d by d matrix of integers giving winding vectors as columns
    # This function should return a full Lattice object

    dim = size(avecs,1)
    int_pts = bravais_int(Xvecs) # integer coordinates of distinct lattice points
    nuc = size(int_pts,2) # number of unit cells
    nbasis = size(basis,2) # number of sites in the primitive unit cell
    nsites = nuc*nbasis # total number of sites in the cluster
   
    # form real space winding vectors
    wind = avecs*Xvecs
    
    # form reciprocal space vectors for supercell
    if dim == 3
    	newwind = wind
    elseif dim == 2
    	newwind = zeros(3,3)
    	newwind[1:2,1:2] = wind
    	newwind[3,3] = 1.0
    else
        error("Lattice must have dimension 2 or 3")
    end
    tp = dot(newwind[:,1],cross(newwind[:,2],newwind[:,3]))
    bigrecip = Array{Float64}(3,3)
    for a in 1:3
        bigrecip[:,a] = cross(newwind[:,mod(a,3)+1],newwind[:,mod(a+1,3)+1])/tp
    end
    if dim == 3
        recip = bigrecip
    elseif dim == 2
        recip = bigrecip[1:2,1:2]
    end
    
    function reduced_pt(pt)
        # gives a unique point inside a supercell defined by wind
        fracpt = mod(recip'*pt,1)
        newpt = wind*fracpt
        return newpt
    end

    function distance(pt1,pt2)
        # finds closest distance between two points taking into account periodicity
        # round to 3 digits
        diff = pt1-pt2
        fracdiff = recip'*diff
        normfrac = mod(fracdiff+0.5*ones(dim),1)-0.5*ones(dim)
        newdiff = wind*normfrac
        return round(norm(newdiff),3)
    end

    # Compute all the sites in the unit cell
    bravais_points = avecs*int_pts
    sites = Array{Float64}(dim,nuc*nbasis) # real space coordinate of site i
    flavors = Array{Int64}(nuc*nbasis) # basis type of site i
    for i in 1:nuc
        for a in 1:nbasis
            sites[:,nuc*(a-1)+i] = reduced_pt(bravais_points[:,i]+basis[:,a]) 
            flavors[nuc*(a-1)+i] = a
        end
    end

    # Compute neighbors
    # For the moment, make simple assumptions:
    # Find only first and second neighbors
    # Assume the same number of these for every site
    # Store them in nns[a,i] and nnns[a,i], here i is the site, and a runs over the nns and nnns of that site

    # First figure out how many there are and at what distance
    alldist = Array{Float64}(0)
    for i in 2:nsites
        push!(alldist,distance(sites[:,1],sites[:,i]))
    end
    dists = unique(alldist)
    sort!(dists)
    nndist = dists[1]
    nnndist = dists[2]
    numnns = length(find(x->x==nndist,alldist))
    numnnns = length(find(x->x==nnndist,alldist))

    # Now construct all of them and store in two 2d arrays
    thenns = Array{Int64}(numnns,nsites)
    thennns = Array{Int64}(numnnns,nsites)
    for i in 1:nsites
        dists = [distance(sites[:,i],sites[:,j]) for j in 1:nsites]
        thenns[:,i] = find(x->x==nndist,dists)
        thennns[:,i] = find(x->x==nnndist,dists)
    end

    # Also compute adjacency matrix: entries = 1 if NN, 2 if NNN, 0 otherwise
    adjacency = Array{Int64}(nsites,nsites)
    for i in 1:nsites
        for j in 1:nsites
            dist = distance(sites[:,i],sites[:,j])
            if dist == nndist
                adjacency[i,j]=1
            elseif dist == nnndist
                adjacency[i,j]=2
            else
                adjacency[i,j]=0
            end
        end
    end

    Lattice(sites,flavors,thenns,thennns,adjacency,nsites,nuc,numnns,numnnns)
end

end
