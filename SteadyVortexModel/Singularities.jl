
module Singularities


export LinearStrengthLineVortex,LinearStrengthLineSource,induced_velocity,SingularityCollection

abstract type Singularity end
abstract type PointSingularity <: Singularity end 
abstract type LineSingularity <: Singularity end 


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
    with strength γ_start and ending at (x_end,y_end) with strength
    γ_end. 
    """
    x_start::A
    y_start::B
    x_end::C
    y_end::D
    σ_start::E
    σ_end::F
end

Base.Broadcast.broadcastable(a::LinearStrengthLineVortex) = Ref(a)
Base.Broadcast.broadcastable(a::LinearStrengthLineSource) = Ref(a)
Base.Broadcast.broadcastable(a::SingularityCollection) = Ref(a)
mutable struct SingularityCollection
    singularities::Vector{<:Singularity}
end 

singularity_length(singularity::T where T<:LineSingularity) = sqrt((singularity.x_end - singularity.x_start)^2 +  (singularity.y_end - singularity.y_start)^2)

# TODO  make a macro strengths(panel::T where T<:LineSingularity) that finds strengths


function induced_velocity_panel_coordinates(xₚ_target,yₚ_target,source::LinearStrengthLineVortex)
    panel_length = singularity_length(source)
    is_on_panel = abs(yₚ_target) <= 1e-8
    r₁ = sqrt(xₚ_target^2 + yₚ_target^2)
    r₂ = sqrt((xₚ_target - panel_length)^2 + yₚ_target^2)
    is_on_endpoint = (r₂ == 0) || (r₁ == 0)
    r₁ = r₁ == 0 ? 1 : r₁
    r₂ = r₂ == 0 ? 1 : r₂
    ln_r₂_r₁ = log(r₂/r₁)
    θ₁ = atan(yₚ_target,xₚ_target)
    θ₂ = atan(yₚ_target,xₚ_target-panel_length)
    dθ = θ₂ - θ₁
    dγ = source.γ_end - source.γ_start
    yₚ_target_regularized = is_on_panel ? 1 : yₚ_target

    uₚ_term1 = yₚ_target / (2π) * dγ / panel_length * ln_r₂_r₁
    uₚ_term2 = (source.γ_start * panel_length + dγ * xₚ_target) / (2π * panel_length) * dθ
    uₚ = is_on_panel ? 0 : uₚ_term1 + uₚ_term2

    vₚ_term1 = (source.γ_start * panel_length + dγ * xₚ_target) / (2π * panel_length) * ln_r₂_r₁
    vₚ_term2 = is_on_panel ? dγ / (2π) : yₚ_target / (2π) * dγ / panel_length * (panel_length / yₚ_target_regularized - dθ)
    vₚ = vₚ_term1 + vₚ_term2

    uₚ = is_on_endpoint ? 0 : uₚ
    vₚ = is_on_endpoint ? 0 : vₚ
    
    return uₚ,vₚ
end 

function induced_velocity_panel_coordinates(xₚ_target,yₚ_target,source::LinearStrengthLineSource)
    panel_length = singularity_length(source)
    is_on_panel = abs(yₚ_target) <= 1e-8
    r₁ = sqrt(xₚ_target^2 + yₚ_target^2)
    r₂ = sqrt((xₚ_target - panel_length)^2 + yₚ_target^2)
    is_on_endpoint = (r₂ == 0) || (r₁ == 0)
    r₁ = r₁ == 0 ? 1 : r₁
    r₂ = r₂ == 0 ? 1 : r₂
    ln_r₂_r₁ = log(r₂/r₁)
    θ₁ = atan(yₚ_target,xₚ_target)
    θ₂ = atan(yₚ_target,xₚ_target-panel_length)
    dθ = θ₂ - θ₁
    dσ = source.σ_end - source.σ_start
    yₚ_target_regularized = is_on_panel ? 1 : yₚ_target

    vₚ_term1 = yₚ_target / (2π) * dσ / panel_length * ln_r₂_r₁
    vₚ_term2 = (source.σ_start * panel_length + dσ * xₚ_target) / (2π * panel_length) * dθ
    vₚ = is_on_panel ? 0 : vₚ_term1 + vₚ_term2

    uₚ_term1 = -(source.σ_start * panel_length + dσ * xₚ_target) / (2π * panel_length) * ln_r₂_r₁
    uₚ_term2 = is_on_panel ? -dσ / (2π) : - yₚ_target / (2π) * dσ / panel_length * (panel_length / yₚ_target_regularized - dθ)
    uₚ = uₚ_term1 + uₚ_term2

    uₚ = is_on_endpoint ? 0 : uₚ
    vₚ = is_on_endpoint ? 0 : vₚ

    return uₚ,vₚ
end 

function induced_velocity(x_target,y_target,source::LinearStrengthLineVortex)
    dx_panel = source.x_end - source.x_start
    dy_panel = source.y_end - source.y_start
    panel_length = sqrt(dx_panel^2 + dy_panel^2)
    xₚ_hat_x = dx_panel/panel_length
    xₚ_hat_y = dy_panel/panel_length
    yₚ_hat_x = - xₚ_hat_y
    yₚ_hat_y = xₚ_hat_x
    x_target_relative = x_target - source.x_start
    y_target_relative = y_target - source.y_start
    xₚ_target = x_target_relative * xₚ_hat_x + y_target_relative * xₚ_hat_y
    yₚ_target = y_target_relative * yₚ_hat_x + y_target_relative * yₚ_hat_y
    uₚ,vₚ = induced_velocity_panel_coordinates(xₚ_target,yₚ_target,source)
    u = uₚ * xₚ_hat_x + vₚ * yₚ_hat_x
    v = uₚ * xₚ_hat_y + vₚ * yₚ_hat_y
    println(u)
    return u,v
end

function induced_velocity(x_target,y_target,source::LinearStrengthLineSource)
    dx_panel = source.x_end - source.x_start
    dy_panel = source.y_end - source.y_start
    panel_length = sqrt(dx_panel^2 + dy_panel^2)
    xₚ_hat_x = dx_panel/panel_length
    xₚ_hat_y = dy_panel/panel_length
    yₚ_hat_x = - xₚ_hat_y
    yₚ_hat_y = xₚ_hat_x
    x_target_relative = x_target - source.x_start
    y_target_relative = y_target - source.y_start
    xₚ_target = x_target_relative * xₚ_hat_x + y_target_relative * xₚ_hat_y
    yₚ_target = y_target_relative * yₚ_hat_x + y_target_relative * yₚ_hat_y
    uₚ,vₚ = induced_velocity_panel_coordinates(xₚ_target,yₚ_target,source)
    u = uₚ * xₚ_hat_x + vₚ * yₚ_hat_x
    v = uₚ * xₚ_hat_y + vₚ * yₚ_hat_y
    return u,v
end

function induced_velocity(x_target,y_target,source::SingularityCollection)
    u,v = 0,0
    for singularity in source.singularities
        u,v = induced_velocity(x_target,y_target,singularity)
    end 
    return u,v
end 



end 





