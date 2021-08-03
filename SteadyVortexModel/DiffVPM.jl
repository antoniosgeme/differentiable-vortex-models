#= 
This code computes the flow over an airfoil with
either a sharp or blunt trailing edge using linearly varying
line vortex panels and a linearly varying source panel if 
needed for the trailing edge. Each panel consists of a line singularity
across it as well as a collocation point in the middle
=#

include("AirfoilModule.jl")
using .AirfoilModule

abstract type Singularity end
abstract type PointSingularity <: Singularity end 
abstract type LineSingularity <: Singularity end 
abstract type Geometry end 
abstract type Point <: Geometry end
abstract type Line <: Geometry end

struct Panel{A,B,C,D,E,F,G,H,I} <: Line
    x_start::A
    y_start::B
    x_end::C
    y_end::D
    x_normal::E
    y_normal::F
    x_tangent::G
    y_tangent::H
    singularity::I
end 

get_midpoint(p::T where T <: Line) = (p.x_end + p.x_start)/2, (p.y_end + p.y_start)/2 

struct LinearStrengthLineVortex{A,B,C,D,E,F} <: LineSingularity
    """
    A linear strength line vortex starting at (x_start,y_start)
    with strength γ_start and ending at (x_end,y_end) with strength
    γ_end. 
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
    A linear strength line vortex starting at (x_start,y_start)
    with strength σ_start and ending at (x_end,y_end) with strength
    σ_end. 
    """
    x_start::A
    y_start::B
    x_end::C
    y_end::D
    σ_start::E
    σ_end::F
end

function induced_velocity(target::Panel,source::LinearStrengthLineVortex; unit_forcing=true)
    dx_panel = source.x_end - source.x_start
    dy_panel = source.y_end - source.y_start
    panel_length = hypot(dx_panel , dy_panel)
    xₚ_hat_x = dx_panel/panel_length
    xₚ_hat_y = dy_panel/panel_length
    yₚ_hat_x = - xₚ_hat_y
    yₚ_hat_y = xₚ_hat_x
    x_target_relative = target.x_target - source.x_start
    y_target_relative = target.y_target - source.y_start
    xₚ_target = x_target_relative * xₚ_hat_x + y_target_relative * xₚ_hat_y
    yₚ_target = y_target_relative * yₚ_hat_x + y_target_relative * yₚ_hat_y

    r₁ = hypot(xₚ_target, yₚ_target)
    r₂ = hypot((xₚ_target - panel_length), yₚ_target)
    ln_r₂_r₁ = log(r₂/r₁)
    θ₁ = atan(yₚ_target,xₚ_target)
    θ₂ = atan(yₚ_target,xₚ_target-panel_length)
    dθ = θ₂ - θ₁

    γ₁,γ₂ = unit_forcing ? (1,1) : (source.γ_start,source.γ_end)

    uₚᵃ = (dθ*(panel_length - xₚ_target) - yₚ_target * ln_r₂_r₁) / (2π * panel_length) * γ₁
    vₚᵃ = ((panel_length - xₚ_target) * ln_r₂_r₁ - panel_length + yₚ_target * dθ) / (2π * panel_length) * γ₁
    uᵃ = uₚᵃ * xₚ_hat_x + vₚᵃ * yₚ_hat_x
    vᵃ = uₚᵃ * xₚ_hat_y + vₚᵃ * yₚ_hat_y
    
    uₚᵇ = (yₚ_target * ln_r₂_r₁  + dθ * xₚ_target) / (2π * panel_length) * γ₂
    vₚᵇ = (xₚ_target * ln_r₂_r₁ + panel_length - yₚ_target * dθ) / (2π * panel_length) * γ₂
    uᵇ = uₚᵇ * xₚ_hat_x + vₚᵇ * yₚ_hat_x
    vᵇ = uₚᵇ * xₚ_hat_y + vₚᵇ * yₚ_hat_y

    return uᵃ,vᵃ,uᵇ,vᵇ
end 


function create_panel(x_start,y_start,x_end,y_end)
