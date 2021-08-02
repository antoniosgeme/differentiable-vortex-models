#= 
This code computes the steady flow over an airfoil with
either a sharp or blunt trailing edge using linearly varying
line vortex panels and a linearly varying source panel if 
needed for the trailing edge. The problem is cast as an optimization
problem and solved using Ipopt. 
=#

include("AirfoilModule.jl")
using .AirfoilModule

abstract type Singularity end
abstract type PointSingularity <: Singularity end 
abstract type LineSingularity <: Singularity end 
abstract type Geometry end 
abstract type Point <: Geometry end
abstract type Line <: Geometry end

struct Panel{A,B,C,D,E,F,G,H,I,J} <: Line
    x_start::A
    y_start::B
    x_end::C
    y_end::D
    x_normal::E
    y_normal::F
    x_tangent::G
    y_tangent::H
    strength_start::I
    strength_end::J
end 

get_midpoint(p::Panel) = (x_end + x_start)/2, (y_end + y_start)/2 

struct LinearStrengthLineVortex{A,B,C,D,E,F} <: LineSingularity
    """ 
    A linear strength line vortex starting at (x_start,y_start)
    with strength γ_start and ending at (x_end,y_end) with strength
    γ_end. It can be attached to any panel and will inherit its location properties.
    """
    x_start::A
    y_start::B
    x_end::C
    y_end::D
    γ_start::E 
    γ_end::F
end  

struct LinearStrengthLineSource{A,B,C,D,E,F} <: LineSingularity
    """
    A linear strength line source starting at (x_start,y_start)
    with strength σ_start and ending at (x_end,y_end) with strength
    σ_end. It can be attached to any panel and will inherit its location properties.
    """
    x_start::A
    y_start::B
    x_end::C
    y_end::D
    
    σ_start::E
    σ_end::F
end

struct CollocationPoint{T} <: Point
    """
    Location where no through flow is enforced on a body. 
    It can be attached to any panel and will inherit its geometric properties.
    """
    x::T
    y::T
    
end 

function create_panel(x_start,y_start,x_end,y_end)
