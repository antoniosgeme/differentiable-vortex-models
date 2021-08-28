include("AirfoilModule.jl")
using .AirfoilModule

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
    dγ = γ₂ - γ₁

    uₚ = yₚ_target / (2π) * dγ / panel_length * ln_r₂_r₁ + (γ₁ * panel_length + dγ * xₚ_target) * dθ / (2π * panel_length) 
    vₚ = (γ₁ * panel_length + dγ * xₚ_target) * ln_r₂_r₁ / (2π * panel_length) + dγ / (2π * panel_length) * (panel_length + yₚ_target - dθ)
    u = uₚ * xₚ_hat_x + vₚ * yₚ_hat_x
    v = uₚ * xₚ_hat_y + vₚ * yₚ_hat_y
    return u,v
end 