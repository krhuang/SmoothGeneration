using Test
# The GH repository implementation https://github.com/uncrayon/Convex-Hull-Jarvis/blob/master/Jarvis.ipynb produces bugs
# Here I implemented it myself, with help from this resource: https://www.geeksforgeeks.org/convex-hull-using-jarvis-algorithm-or-wrapping/
function orientation(p, q, r)
    value = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])
    if value == 0
        return 0 # Colinear
    elseif value < 0
        return -1 # Counter-clockwise
    else 
        return 1 # Clockwise
    end
end

# Convex_Hull takes as input a nx2 array of integers, representing points in the plane. It returns the vertices of the convex hull, in clockwise position. 
function Convex_Hull(points)
    output_array = Vector{Vector{Int}}()
    # Initialize the left-most point
    starting_point = points[1, :]
    for i in 2:size(points, 1)
        if points[i, 1] < starting_point[1] || 
            (points[i, 1] == starting_point[1] && points[i, 2] < starting_point[2])
            starting_point = points[i, :]
        end
    end
    push!(output_array, starting_point)
    current_point = starting_point
    while(true)
        next_point = points[1, :] # Grab something as the next point. If it happens to be the same as the current point, then just grab a different one.
        if next_point == current_point
            next_point = points[2, :]
        end
        for point in eachrow(points)
            if orientation(current_point, next_point, point) == -1 # Counter-clockwise orientation or co-linear
                next_point = point
            end
        end
        if next_point == starting_point
            break
        end
        push!(output_array, next_point)
        
        current_point = next_point
    end
    return output_array
end
     

# TestCase
points = [0 3; 2 2; 1 1; 2 1; 3 0; 0 0; 3 3]
hull = Convex_Hull(points)
println(hull)