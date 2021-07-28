module AirfoilModule

export Airfoil,get_upper_coordinates,get_lower_coordinates,scale!,get_area,get_chord,get_centroid,set_angle_of_attack!,write_file,repanel!,sharp_trailing_edge

using Dierckx
using Printf

mutable struct Airfoil
    "Name of the airfoil eg naca0012"
    name::String
    "x coordinates of the airfoil"
    x::Array{Float64,1}
    "y coordinates of the airfoil"
    y::Array{Float64,1}
    "angle of attack of the airfoil"
    angle_of_attack::Real
end 

Base.show(io::IO,airfoil::Airfoil) = println(io,"I am a "*airfoil.name*" airfoil, at an angle of attack of "*string(airfoil.angle_of_attack)*"°!")


function Airfoil(name::String; points_per_side =50, filepath::String = "",sharp_trailing_edge=false,fine_trailing_edge_paneling=true)
    if isempty(filepath)
        if occursin("naca",lowercase(name))
            try
                x,y = naca4(name,sharp_trailing_edge=sharp_trailing_edge,fine_trailing_edge_paneling=fine_trailing_edge_paneling)
                airfoil = Airfoil(name,x,y,0)
            catch 
                x,y = UIUC_coords(name)
                airfoil = Airfoil(name,x,y,0)
            end 
        else
            x,y = UIUC_coords(name)
            airfoil = Airfoil(name,x,y,0)
        end
    else
        x,y = from_filepath(filepath)
        airfoil = Airfoil(name,x,y,0)
    end

    return airfoil
end


"""
    naca4(name::String, points_per_side::Int64 = 50;fine_trailing_edge_paneling=false)

Generates coordinates for 4 digit NACA airfoils with 2*points_per_side-1 points.
If fine_trailing_edge_paneling=false cosine spacing is only used on the LE while fine_trailing_edge_paneling=true also
applies cosine spacing to the TE.
"""
function naca4(name::String, points_per_side::Int64 = 50; sharp_trailing_edge=false,fine_trailing_edge_paneling=false)

    naca_num = strip(name)[5:end]

    if length(naca_num) != 4
        Error("Oops! Can only populate from 4 digit naca airfoils")
    end

    m = parse(Float64,string(naca_num[1]))/100.
    p = parse(Float64,string(naca_num[2]))/10.
    t = parse(Float64,string(naca_num[3:end]))/100.

    fine_trailing_edge_paneling ? x_t = cos_space(0,1,points_per_side) : x_t = half_cos_space(points_per_side)[end:-1:1]

    x_t1 = x_t[x_t .<= p]
    x_t2 = x_t[x_t .> p]

    TE_coefficient = sharp_trailing_edge ? 0.1036 : 0.1015 
    y_t(x) = 5t*(0.2969 * √x - 0.1260x - 0.3516x^2 + 0.2843x^3 - TE_coefficient*x^4) # 0.1015/0.1036 for blunt/sharp TE

    p == 0  ? p = 0.5 : nothing

    y_c1(x) =(m / p^2) * (2p * (x) - x^2)
    y_c2(x) = m / (1 - p)^2 * ((1 - 2p) + 2p * x - x^2)

    y_c = vcat(y_c1.(x_t1) , y_c1.(x_t2))

    dyc1_dx(x) = 2m / p^2 * (p-x)
    dyc2_dx(x) = 2m / (1-p)^2 * (p-x)

    dyc_dc = vcat(dyc1_dx.(x_t1),dyc2_dx.(x_t2))

    θ = atan.(dyc_dc)

    x_U = similar(x_t)
    x_L = similar(x_t)
    y_U = similar(x_t)
    y_L = similar(x_t)

    @. x_U = x_t - y_t(x_t) * sin(θ)
    @. x_L = x_t + y_t(x_t) * sin(θ)
    @. y_U = y_c + y_t(x_t) * cos(θ)
    @. y_L = y_c - y_t(x_t) * cos(θ)

    # Flip upper surface so it's back to front
    x_U = x_U[end:-1:1, :]
    y_U = y_U[end:-1:1, :]

    # Trim 1 point from lower surface so there's no overlap
    x_L = x_L[2:end]
    y_L = y_L[2:end]

    x = vcat(x_U, x_L)
    y = vcat(y_U, y_L)

    return vec(x),vec(y)
end

"""
Internal function used to populate airfoil coordinates from UIUC database.
if UIUC database is not in current directory specify the location of the UIUC folder
through the argument filepath. 
"""
function UIUC_coords(name,filepath::String="")

    if isempty(filepath)
        datfile = strip(lowercase(name))*".dat"
        dir = raw".\\airfoil_database"
        
    else
        dir = filepath*raw".\\airfoil_database"
    end 

    cd(dir)

    if datfile in readdir()
         println("Airfoil found!")

         f = open(datfile,"r")
         io = readlines(f)
         close(f)
         cd(raw"..\\..\\")

         x = [parse(Float64,split(xx)[1]) for xx in io[2:end]]
         y = [parse(Float64,split(yy)[2]) for yy in io[2:end]]

         return x,y

    else
        cd(raw"..\\..\\")
        println("No such airfoil exists in the UIUC database")
        return NaN,NaN
    end
end

"""
Internal function used to populate airfoil coordinates from fileapath.
"""
function from_filepath(name)

    if  ! occursin(".dat",name)
        datfile = strip(name)*".dat"
    else
        datfile = strip(name)
    end

    f = open(datfile,"r")

    io = readlines(f)

    close(f)

    x = [parse(Float64,split(xx)[1]) for xx in io[2:end]]
    y = [parse(Float64,split(yy)[2]) for yy in io[2:end]]

    return x,y
end

function  get_upper_coordinates(airfoil::Airfoil) 
    x = airfoil.x[argmin(airfoil.x):-1:1]
    y = airfoil.y[argmin(airfoil.x):-1:1]
    return x,y
end 

function get_lower_coordinates(airfoil::Airfoil) 
    x = airfoil.x[argmin(airfoil.x[:,1]):end]
    y = airfoil.y[argmin(airfoil.x[:,1]):end]
    return x,y
end 


function scale!(airfoil::Airfoil,scalefactor::Real) 
    airfoil.x .*= scalefactor  
    airfoil.y .*= scalefactor
end 

get_area(airfoil::Airfoil) = 0.5 * sum(airfoil.x .* circshift(airfoil.y,-1)
                                 .- airfoil.y .* circshift(airfoil.x,-1))


function get_chord(airfoil::Airfoil)
    temp = deepcopy(airfoil)
    setAoA!(temp,0)
    chord = maximum(temp.x) - minimum(temp.x)
    return chord
end

function get_centroid(airfoil::Airfoil)
    x =  airfoil.x
    y  = airfoil.y

    xn = circshift(airfoil.x,-1)
    yn = circshift(airfoil.y,-1)

    a = x .* yn .- y .* xn

    A = 0.5 * sum(a)

    x_c = 1 / (6 * A) * sum(a .* (x .+ xn))
    y_c = 1 / (6 * A) * sum(a .* (y .+ yn))

    return x_c , y_c
end


function repanel!(airfoil::Airfoil,points_per_side;fine_trailing_edge_paneling=false)

    if occursin("naca",lowercase(airfoil.name)) && length(airfoil.name[5:end]) == 4 
        airfoil.x, airfoil.y = naca4(airfoil.name,points_per_side,fine_trailing_edge_paneling=fine_trailing_edge_paneling)
        if airfoil.angle_of_attack != 0
            target_aoa = airfoil.angle_of_attack
            airfoil.angle_of_attack = 0
            set_angle_of_attack(airfoil,target_aoa)
        end

    else
        upper_x, upper_y = get_upper_coordinates(airfoil)
        lower_x, lower_y = get_lower_coordinates(airfoil)

        interp_upper = Spline1D(upper_x,upper_y)

        interp_lower = Spline1D(lower_x,lower_y)

        if fine_trailing_edge_paneling
            s = cos_space(0,1,points_per_side)
        else
            s = half_cos_space(points_per_side)
        end

        new_upper = interp_upper.(s)
        new_lower = interp_lower.(s)

        new_x = vcat(s[end:-1:1],s[2:end])

        new_y = vcat(new_upper[end:-1:1],new_lower[2:end])

        airfoil.x,airfoil.y = new_x,new_y
    end
end

function set_angle_of_attack!(airfoil::Airfoil,angle_of_attack::Real)
    AoA_diff = angle_of_attack - airfoil.angle_of_attack
    AoA_diff = -deg2rad(AoA_diff)
    for i in 1:length(airfoil.x)
         airfoil.x[i],airfoil.y[i] = rotate2D([airfoil.x[i],airfoil.y[i]],AoA_diff)
    end
    airfoil.angle_of_attack = angle_of_attack
end

function write_file(airfoil::Airfoil)
    open(airfoil.name*".dat","w") do f
        write(f,airfoil.name*"\n")
        for i in 1:length(airfoil.x)
            var = (airfoil.x[i],airfoil.y[i])
            str = "$(@sprintf("     %.12f    %.12f\n",var...))"
            write(f,str)
        end
    end
end

sharp_trailing_edge(airfoil::Airfoil) =  isapprox(0,airfoil.y[1]-airfoil.y[end],atol=1e-10)


cos_space(min,max,npoints) = (max+min)/2 .+ (max-min)/2 * cos.(Array{Float64}(LinRange(π, 0, npoints)))

half_cos_space(npoints) = 1 .- cos.(LinRange(π, 0, npoints) ./ 2)

"Rotates 2D vector v by angle θ (must be in radians)"
rotate2D(v,θ) = [ cos(θ) -sin(θ)
                  sin(θ)  cos(θ) ] * v


end 
