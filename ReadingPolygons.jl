using HDF5
using Polymake 	#via polymake.jl
using Test 		#for test cases
include("jarvismarch.jl")


function read_numbers_of_polygons(file::String)
    # Open the HDF5 file in read-only mode
    numbers_of_polygons = []
    h5open(file, "r") do f
        # Read the "numbers_of_polygons" dataset
        numbers_of_polygons = read(f["numbers_of_polygons"])
    end
    return numbers_of_polygons
end

function read_polygons(file::String, dataset_name::String)
    polygons = []
    h5open(file, "r") do f
        if haskey(f, dataset_name)
            # Read the dataset containing polygons
            polygons = read(f[dataset_name])
        else
            println("Dataset $dataset_name not found in file $file")
        end
    end
    return polygons
end

function print_polygons(polygons, n_vertices::Int) # Removing type requirements on polygons fixed it, for some reason
    # Each polygon is stored as a NamedTuple with 2n fields
    println("Number of polygons with $n_vertices vertices: $(length(polygons))")
    for (i, polygon) in enumerate(polygons)
        println("Polygon $i:")
        for j in 1:n_vertices
            x = getfield(polygon, Symbol(string(2 * j - 1)))
            y = getfield(polygon, Symbol(string(2 * j)))
            print(" ($x, $y)")
        end
        println()
    end
end

function process_polygons(polygons, n_vertices::Int) 	# Removing type requirements on polygons fixed it, for some reason
														# Similar to print_polygons but instead checks if they are smooth via polymake, and counts their #interior lattice points
    # Each polygon is stored as a NamedTuple with 2n fields
    println("Number of polygons with $n_vertices vertices: $(length(polygons))")

    for (i, polygon) in enumerate(polygons)
    	vertices_matrix = zeros(Int, n_vertices, 2)
    	polymake_vertices_matrix = zeros(Int, n_vertices, 3)
        for j in 1:n_vertices
            x = getfield(polygon, Symbol(string(2 * j - 1)))
            y = getfield(polygon, Symbol(string(2 * j)))
            vertex_vector = [x; y]
            vertices_matrix[j, :] .= vertex_vector

            polymake_vertex_vector = [1; x; y] # Appending 1 due to Polymake standards
            polymake_vertices_matrix[j, :] = polymake_vertex_vector
            
        end
        new_polygon = polytope.Polytope(POINTS=polymake_vertices_matrix)
        if new_polygon.SMOOTH
        	print(n_vertices, ":")
        	print(new_polygon.N_INTERIOR_LATTICE_POINTS, " ")
        	println(join(vec(vertices_matrix), " "))
        end
    end
end

function process_hdf5_file(file::String)
    println("Processing file: $file")
    
    # Step 1: Read the "numbers_of_polygons" dataset
    numbers_of_polygons = read_numbers_of_polygons(file)
    println("Numbers of polygons per vertex count:")
    println(numbers_of_polygons)

    # Step 2: Iterate through the datasets n3, n4, ..., and read polygons
    for (i, count) in enumerate(numbers_of_polygons)
        n_vertices = i + 2  # Dataset starts at n3
        dataset_name = "n$n_vertices"
        println("\nReading dataset: $dataset_name")
        
        polygons = read_polygons(file, dataset_name)
        if length(polygons) > 0
            process_polygons(polygons, n_vertices)
        else
            println("No polygons found in dataset $dataset_name.")
        end
    end
end

# Main script
file_name = "koelman70/l16.h5"  # Replace with the actual file name
process_hdf5_file(file_name)