include("AirfoilModule.jl")
using .AirfoilModule

#using Plots
#plotly()
import PyPlot as plt

abstract type Point end
abstract type Singularity end 

struct CollocationPoint{T} <: Point
    x::T
    y::T
    x_normal::T
    y_normal::T
    x_tangent::T
    y_tangent::T
end 


struct VortexPoint{T} <: Singularity
    x::T
    y::T
    x_normal::T
    y_normal::T
    x_tangent::T
    y_tangent::T
    γ::T
end 

struct SourcePoint{T} <: Singularity
    x::T
    y::T
    x_normal::T
    y_normal::T
    x_tangent::T
    y_tangent::T
    σ::T
end 


struct Body 
    singularities::Array{<: Singularity,1}
    collocations::Array{CollocationPoint,1}
    num_panels::Int
end 

struct OperatingPoint 
    U_freestream::Float64
    V_freestream::Float64
end 

OperatingPoint(u) = OperatingPoint(u,0)
VortexPoint(a,b,c,d,e,f) = VortexPoint(a,b,c,d,e,f,1)
SourcePoint(a,b,c,d,e,f) = SourcePoint(a,b,c,d,e,f,1)
 
function induced_velocity(target::CollocationPoint,source::VortexPoint; unit_forcing=true)
        A = [0 1
            -1 0]
    X = [target.x - source.x ,target.y - source.y ]
    r² = (target.x - source.x)^2 + (target.y - source.y)^2
    γ = unit_forcing ? 1 : source.γ
    C = γ/(2π*r²)
    u,w = C*A*X
    return u,w
end 

function induced_velocity(target::CollocationPoint,source::SourcePoint; unit_forcing=true)
        A = [0 1
            1 0]
    X = [target.y - source.y ,target.x - source.x ]
    r² = (target.x - source.x)^2 + (target.y - source.y)^2
    σ = unit_forcing ? 1 : source.σ
    C = σ/(2π*r²)
    u,w = C*A*X
    return u,w
end 

function populate_influence_matrix(body::Body)
    A = zeros(body.num_panels+1,body.num_panels+1)
    for (i,collocation) in enumerate(body.collocations)
        for (j,singularity) in enumerate(body.singularities)
            u,w = induced_velocity(collocation,singularity)
            A[i,j] = dot([u,w],[collocation.x_normal,collocation.y_normal])
        end 
    end 
    A[end,1] = 1
    A[end,end] = 1
    return A
end 

function populate_RHS(body::Body,op::OperatingPoint)
    U_freestream = op.U_freestream
    V_freestream = op.V_freestream
    RHS = zeros(body.num_panels + 1)
    for (i,collocation) in enumerate(body.collocations)
        RHS[i] = -dot([U_freestream, V_freestream],[collocation.x_normal,collocation.y_normal])
    end 
    return RHS
end

function enforce_no_throughflow!(body::Body,op::OperatingPoint)
    A = populate_influence_matrix(body)
    RHS = populate_RHS(body,op)
    
    Γ = A \ RHS
    return Γ
end 

function construct_geometry(airfoil::Airfoil)
    dx_panels = airfoil.x - circshift(airfoil.x,-1)
    dy_panels = airfoil.y - circshift(airfoil.y,-1)
    sharp_trailing_edge(airfoil) ? begin pop!(dx_panels)
                                         pop!(dy_panels) 
                                    end : nothing     
    
    panel_lengths = hypot.(dx_panels , dy_panels)
    x_normal = -dy_panels ./ panel_lengths
    y_normal = dx_panels ./ panel_lengths
    y_tangent = - copy(x_normal)
    x_tangent = copy(y_normal) 
    x_midpoints, y_midpoints = compute_midpoints(airfoil)
    collocations = CollocationPoint.(   
                                    x_midpoints,
                                    y_midpoints,
                                    x_normal,
                                    y_normal,
                                    x_tangent,
                                    y_tangent
                                    )
    if sharp_trailing_edge(airfoil)
        vortices = VortexPoint.(
                                airfoil.x,
                                airfoil.y,
                                push!(x_normal,x_normal[end]),
                                push!(y_normal,y_normal[end]),
                                push!(x_tangent,x_tangent[end]),
                                push!(y_tangent,y_tangent[end]),
                                )  
    else 
        vortices = VortexPoint.(airfoil.x,
                                airfoil.y,
                                push!(x_normal[1:end-1],x_normal[end-1]),
                                push!(y_normal[1:end-1],y_normal[end-1]),
                                push!(x_tangent[1:end-1],x_tangent[end-1]),
                                push!(y_tangent[1:end-1],y_tangent[end-1]),
                                )  
        #= 
        source = SourcePoint(
                            x_midpoints[end],
                            y_midpoints[end],
                            x_normal[end],
                            y_normal[end],
                            x_tangent[end],
                            y_tangent[end])
        =#
        
    end 


    return x_normal, y_normal, x_tangent, y_tangent
end 


function compute_midpoints(airfoil::Airfoil)
    if sharp_trailing_edge(airfoil)
        x_midpoints = (airfoil.x[2:end] +  airfoil.x[1:end-1])/2  
        y_midpoints = (airfoil.y[2:end] +  airfoil.y[1:end-1])/2
    else 
        x_midpoints = (airfoil.x + circshift(airfoil.x,-1))/2 
        y_midpoints = (airfoil.y + circshift(airfoil.y,-1))/2 
    end 
    return x_midpoints, y_midpoints
end 






naca = Airfoil("naca6409",sharp_trailing_edge=true)

x_normal, y_normal, x_tangent, y_tangent = construct_geometry(naca)
x_midpoints, y_midpoints = compute_collocation_points(naca)

##

plt.figure()
plt.plot(naca.x,naca.y)
plt.quiver(x_midpoints,y_midpoints,x_normal,y_normal,scale=20)
#plt.quiver(x_midpoints,y_midpoints,x_tangent,y_tangent,scale=20)
plt.axis("equal")

##