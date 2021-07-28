using Plots
plotly()

function Plots.plot(airfoil::Airfoil)
    plot(airfoil.x,airfoil.y,aspect_ratio=:equal,legend=false)
end 

function Plots.scatter(airfoil::Airfoil,ms=3)
    scatter(airfoil.x,airfoil.y,aspect_ratio=:equal,legend=false,ms=ms)
end 


using .Singularities
using .AirfoilModule


airfoil = Airfoil("naca6409")
repanel!(airfoil,10)
set_angle_of_attack!(airfoil,20)
scatter(airfoil)

U = 1 

function compute_collocation_points(airfoil::Airfoil)
    x_midpoints = (airfoil.x + circshift(airfoil.x,-1))/2 
    y_midpoints = (airfoil.y + circshift(airfoil.y,-1))/2 
    return x_midpoints, y_midpoints
end 

x_midpoints,y_midpoints = compute_collocation_points(airfoil)
