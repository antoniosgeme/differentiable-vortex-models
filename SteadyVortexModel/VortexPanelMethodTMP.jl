include("AirfoilModule.jl")
using .AirfoilModule

using Plots
plotly() 

function compute_collocation_points(airfoil::Airfoil)
    x_midpoints = (airfoil.x + circshift(airfoil.x,-1))/2 
    y_midpoints = (airfoil.y + circshift(airfoil.y,-1))/2 
    return x_midpoints, y_midpoints
end 


function VOR2DL(γ₁,γ₂,x_target,y_target,x_start,y_start,x_end,y_end)
    dx_panel = x_end - x_start
    dy_panel = y_end - y_start
    panel_length = hypot(dx_panel , dy_panel)
    xₚ_hat_x = dx_panel/panel_length
    xₚ_hat_y = dy_panel/panel_length
    yₚ_hat_x = - xₚ_hat_y
    yₚ_hat_y = xₚ_hat_x
    x_target_relative = x_target - x_start
    y_target_relative = y_target - y_start
    xₚ_target = x_target_relative * xₚ_hat_x + y_target_relative * xₚ_hat_y
    yₚ_target = y_target_relative * yₚ_hat_x + y_target_relative * yₚ_hat_y

    r₁ = hypot(xₚ_target, yₚ_target)
    r₂ = hypot((xₚ_target - panel_length), yₚ_target)
    ln_r₂_r₁ = log(r₂/r₁)
    θ₁ = atan(yₚ_target,xₚ_target)
    θ₂ = atan(yₚ_target,xₚ_target-panel_length)
    dθ = θ₂ - θ₁

    uₚᵃ = (dθ*(panel_length - xₚ_target) - yₚ_target * ln_r₂_r₁) / (2π * panel_length) * γ₁
    uₚᵇ = (yₚ_target * ln_r₂_r₁  + dθ * xₚ_target) / (2π * panel_length) * γ₂

    vₚᵃ = ((panel_length - xₚ_target) * ln_r₂_r₁ - panel_length + yₚ_target * dθ) / (2π * panel_length) * γ₁
    vₚᵇ = (xₚ_target * ln_r₂_r₁ + panel_length - yₚ_target * dθ) / (2π * panel_length) * γ₂

    uᵃ = uₚᵃ * xₚ_hat_x + vₚᵃ * yₚ_hat_x
    uᵇ = uₚᵇ * xₚ_hat_x + vₚᵇ * yₚ_hat_x

    vᵃ = uₚᵃ * xₚ_hat_y + vₚᵃ * yₚ_hat_y
    vᵇ = uₚᵇ * xₚ_hat_y + vₚᵇ * yₚ_hat_y

    u = uᵃ + uᵇ
    v = vᵃ + vᵇ

    return u,v,uᵃ,vᵃ,uᵇ,vᵇ
end 

function compute_panel_unit_vectors(airfoil::Airfoil)
    dx_panel = diff(airfoil.x)
    dy_panel = diff(airfoil.y)
    panel_length = hypot.(dx_panel , dy_panel)
    xₚ_hat_x = dx_panel ./ panel_length
    xₚ_hat_y = dy_panel ./ panel_length
    yₚ_hat_x = -xₚ_hat_y
    yₚ_hat_y = xₚ_hat_x
    return xₚ_hat_x,xₚ_hat_y,yₚ_hat_x,yₚ_hat_y
end 



airfoil = Airfoil("naca6409")
set_angle_of_attack!(airfoil,20)
scatter(airfoil)
x_midpoints,y_midpoints = compute_collocation_points(airfoil)
xₚ_hat_x,xₚ_hat_y,yₚ_hat_x,yₚ_hat_y = compute_panel_unit_vectors(airfoil)


N = length(airfoil.x)
A = zeros(N+1,N+1)



