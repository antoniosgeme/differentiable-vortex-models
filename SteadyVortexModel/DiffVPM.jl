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
    x_target,y_target = get_midpoint(target)
    dx_panel = source.x_end - source.x_start
    dy_panel = source.y_end - source.y_start
    panel_length = hypot(dx_panel , dy_panel)
    xₚ_hat_x = dx_panel/panel_length
    xₚ_hat_y = dy_panel/panel_length
    yₚ_hat_x = - xₚ_hat_y
    yₚ_hat_y = xₚ_hat_x
    x_target_relative = x_target - source.x_start
    y_target_relative = y_target - source.y_start
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

function induced_velocity(target::Panel,source::LinearStrengthLineSource; unit_forcing=true)
    x_target,y_target = get_midpoint(target)
    dx_panel = source.x_end - source.x_start
    dy_panel = source.y_end - source.y_start
    panel_length = hypot(dx_panel , dy_panel)
    xₚ_hat_x = dx_panel/panel_length
    xₚ_hat_y = dy_panel/panel_length
    yₚ_hat_x = - xₚ_hat_y
    yₚ_hat_y = xₚ_hat_x
    x_target_relative = x_target - source.x_start
    y_target_relative = y_target - source.y_start
    xₚ_target = x_target_relative * xₚ_hat_x + y_target_relative * xₚ_hat_y
    yₚ_target = y_target_relative * yₚ_hat_x + y_target_relative * yₚ_hat_y

    r₁ = hypot(xₚ_target, yₚ_target)
    r₂ = hypot((xₚ_target - panel_length), yₚ_target)
    ln_r₂_r₁ = log(r₂/r₁)
    θ₁ = atan(yₚ_target,xₚ_target)
    θ₂ = atan(yₚ_target,xₚ_target-panel_length)
    dθ = θ₂ - θ₁

    σ₁,σ₂ = unit_forcing ? (1,1) : (source.σ_start,source.σ_end)

    uₚᵃ = - σ₁ *  (x_end - xₚ_target) / (2π * panel_length) * ln_r₂_r₁ +  σ₁ /  (2π * panel_length) * (panel_length  + yₚ_target * dθ)
    vₚᵃ = yₚ_target * σ₁ / (2π * panel_length) * ln_r₂_r₁ + σ₁ *  (x_end - xₚ_target) / (2π * panel_length)  * dθ

    uᵃ = uₚᵃ * xₚ_hat_x + vₚᵃ * yₚ_hat_x
    vᵃ = uₚᵃ * xₚ_hat_y + vₚᵃ * yₚ_hat_y

    uₚᵇ = - σ₂ *  (xₚ_target - x_start) / (2π * panel_length) * ln_r₂_r₁ -  σ₂ /  (2π * panel_length) * (panel_length  + yₚ_target * dθ)
    vₚᵇ = - yₚ_target * σ₂ / (2π * panel_length) * ln_r₂_r₁ + σ₂ *  (xₚ_target - x_start) / (2π * panel_length)  * dθ

    uᵇ = uₚᵇ * xₚ_hat_x + vₚᵇ * yₚ_hat_x
    vᵇ = uₚᵇ * xₚ_hat_y + vₚᵇ * yₚ_hat_y

    return uᵃ,vᵃ,uᵇ,vᵇ
end 

function induced_velocity(target::Panel,source::Panel; unit_forcing=true)
    return induced_velocity(target,source.singularity,unit_forcing=unit_forcing)
end 

function induced_velocity(x_target,y_target,source::Panel; unit_forcing=false,return_total_velocity=true)
    dx_panel = source.x_end - source.x_start
    dy_panel = source.y_end - source.y_start
    panel_length = hypot(dx_panel , dy_panel)
    xₚ_hat_x = dx_panel/panel_length
    xₚ_hat_y = dy_panel/panel_length
    yₚ_hat_x = - xₚ_hat_y
    yₚ_hat_y = xₚ_hat_x
    x_target_relative = x_target - source.x_start
    y_target_relative = y_target - source.y_start
    xₚ_target = x_target_relative * xₚ_hat_x + y_target_relative * xₚ_hat_y
    yₚ_target = y_target_relative * yₚ_hat_x + y_target_relative * yₚ_hat_y

    r₁ = hypot(xₚ_target, yₚ_target)
    r₂ = hypot((xₚ_target - panel_length), yₚ_target)
    ln_r₂_r₁ = log(r₂/r₁)
    θ₁ = atan(yₚ_target,xₚ_target)
    θ₂ = atan(yₚ_target,xₚ_target-panel_length)
    dθ = θ₂ - θ₁

    σ₁,σ₂ = unit_forcing ? (1,1) : (source.σ_start,source.σ_end)

    uₚᵃ = - σ₁ *  (x_end - xₚ_target) / (2π * panel_length) * ln_r₂_r₁ +  σ₁ /  (2π * panel_length) * (panel_length  + yₚ_target * dθ)
    vₚᵃ = yₚ_target * σ₁ / (2π * panel_length) * ln_r₂_r₁ + σ₁ *  (x_end - xₚ_target) / (2π * panel_length)  * dθ

    uᵃ = uₚᵃ * xₚ_hat_x + vₚᵃ * yₚ_hat_x
    vᵃ = uₚᵃ * xₚ_hat_y + vₚᵃ * yₚ_hat_y

    uₚᵇ = - σ₂ *  (xₚ_target - x_start) / (2π * panel_length) * ln_r₂_r₁ -  σ₂ /  (2π * panel_length) * (panel_length  + yₚ_target * dθ)
    vₚᵇ = - yₚ_target * σ₂ / (2π * panel_length) * ln_r₂_r₁ + σ₂ *  (xₚ_target - x_start) / (2π * panel_length)  * dθ

    uᵇ = uₚᵇ * xₚ_hat_x + vₚᵇ * yₚ_hat_x
    vᵇ = uₚᵇ * xₚ_hat_y + vₚᵇ * yₚ_hat_y


    if return_total_velocity
        return uᵃ + uᵇ, vᵃ + vᵇ
    else
        return uᵃ,vᵃ,uᵇ,vᵇ
    end 
end 


function create_panel(x_start,y_start,x_end,y_end;line_vortex=true)
    singularity_type = line_vortex ? LinearStrengthLineVortex : LinearStrengthLineSource
    singularity = singularity_type(x_start,y_start,x_end,y_end,1.,1.)
    dx = x_end - x_start
    dy = y_end - y_start
    panel_length = hypot(dx , dy)
    x_normal = -dy / panel_length
    y_normal = dx / panel_length
    y_tangent = - copy(x_normal)
    x_tangent = copy(y_normal) 
    return Panel(x_start,y_start,x_end,y_end,x_normal, y_normal,x_tangent,y_tangent,singularity)
end

airfoil = Airfoil("naca6409")
panels = create_panel.(airfoil.x[1:end-1],airfoil.y[1:end-1],airfoil.x[2:end],airfoil.x[2:end])