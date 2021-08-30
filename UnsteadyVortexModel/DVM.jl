using Plots
using LinearAlgebra

abstract type Point end 
abstract type Body end  

struct Kinematics
    X₀::Vector{Float64} # Location of pitch axis in inertial frame
    Y₀::Vector{Float64} # Location of pitch axis in inertial frame
    U₀::Vector{Float64} # Velocity of pitch axis in inertial frame
    V₀::Vector{Float64} # Velocity of pitch axis in inertial frame
    α::Vector{Float64} # Vector of the angles of attack
    α̇::Vector{Float64}# Vector of the pitch rate 
end 




mutable struct CollocationPoint{T} <: Point
    x_wing::T
    y_wing::T
    x_inertial::T
    y_inertial::T
    x_inertial_history::Vector{T}
    y_inertial_history::Vector{T}
    x_normal_wing::T
    y_normal_wing::T
    x_normal_inertial::T
    y_normal_inertial::T
    x_tangent_wing::T
    y_tangent_wing::T
    x_tangent_inertial::T
    y_tangent_inertial::T
end 

mutable struct Vortex{T} <: Point
    x_wing::T
    y_wing::T
    x_inertial::T
    y_inertial::T
    x_inertial_history::Vector{T}
    y_inertial_history::Vector{T}
    x_normal_wing::T
    y_normal_wing::T
    x_normal_inertial::T
    y_normal_inertial::T
    x_tangent_wing::T
    y_tangent_wing::T
    x_tangent_inertial::T
    y_tangent_inertial::T
    Γ::T
end 

mutable struct Airfoil <: Body 
    name::String
    chord::Float64
    pitch_axis::Float64 # x coordinate of litching axis in wing frame as a fraction of chord
    x_pitch_axis::Float64 # Location of pitch axis in inertial frame
    y_pitch_axis::Float64 # Location of pitch axis in inertial frame
    vortices::Vector{Vortex}
    collocations::Vector{CollocationPoint}
    leading_edge_shedding::Bool
end 

struct Gust
    u::Vector{Float64}
    v::Vector{Float64}
end 

mutable struct Flowfield
    wing::Airfoil
    wake_vortices::Vector{Vortex}
    gust::Gust
end 

function top_hat(x::Vector{<:Number})
    u = [x[i] > 10 || x[i] < 5 ? 0 : 1 for i in 1:length(x)]
    v = zeros(length(x))
    return Gust(u,v)
end 

function set_kinematics(time::Vector{Float64},U∞::Number,V∞::Number,α::Vector{Float64},airfoil::Airfoil)
    X₀  = U∞.*time .+ airfoil.x_pitch_axis
    Y₀ = V∞.*time .+ airfoil.y_pitch_axis
    U∞ = ones(length(time)) * U∞
    V∞ = ones(length(time)) * V∞
    α̇ = gradient1D(α,time)
    return Kinematics(X₀,Y₀,U∞,V∞,α,α̇ )
end

function gradient1D(f::Vector{<:Number},x::Vector{<:Number})
    df_dx = similar(f)
    for i in 2:length(f)-1
        df_dx[i] = (f[i-1] - f[i+1])/(x[i-1] - x[i+1])
    end 
    df_dx[1] = (f[1] - f[2])/(x[1] - x[2])
    df_dx[end] = (f[end-1] - f[end])/(x[end-1] - x[end])
    return df_dx
end 


"""
    nacacamberline(naca_airfoil::String;chord=1,n_panels=50)

Returns coordinates of the specified 4 digit NACA Airfoil with a 
user defined chord and discretization size. The chord is defaulted to 1 
and the number of panels to 50
# Examples
```julia-repl
julia> x,z = nacacamberline("6400",chord=3,n_panels=100)
julia> x,z = nacacamberline("2600")

```
"""
function nacacamberline(naca_airfoil::String;chord=1,n_panels=50)
    m = parse(Int64,naca_airfoil[1])/100
    p = parse(Int64,naca_airfoil[2])/10
    c = chord
    x = collect(LinRange(0,chord,n_panels))
    y = similar(x)

    if m != 0 && p != 0
        for i in 1:n_panels
            if x[i] <= p*c
                y[i] =  m * (x[i] / (p^2)) * (2p - x[i]/c)
            else 
                y[i] = m * ((c-x[i]) / (1-p)^2) * (1 + x[i] / c - 2p)
            end 
        end 
    end 
    return x,y
end 

function wing_to_inertial_position(x_wing,y_wing,X₀,Y₀,α)
    R = [cosd(α) sind(α)
        -sind(α) cosd(α)]
    
    x_inertial,y_inertial = R*[x_wing,y_wing] + [X₀,Y₀]
    return x_inertial,y_inertial
end 

function wing_to_inertial_velocity(x_wing,y_wing,U∞,W∞,α,α̇ )
    Ṙ = α̇  * [-sind(α)  cosd(α)
              -cosd(α) -sind(α)]
    
    u_inertial,w_inertial = Ṙ * [x_wing,y_wing] + [U∞,W∞]
    return u_inertial,v_inertial
end 

function inertial_to_wing_position(x_inertial,y_inertial,X₀,Y₀,α)
    R = [cosd(α) -sind(α)
         sind(α)  cosd(α)]
    
    x_wing,y_wing = R * [x_inertial - X₀,y_inertial - Y₀] 
    return  x_wing,y_wing
end 

function inertial_to_wing_velocity(x_inertial,y_inertial,X₀,Y₀,U∞,W∞,α,α̇ )
    R = [cosd(α) -sind(α)
         sind(α)  cosd(α)]
    
    Ṙ = α̇  * [-sind(α) -cosd(α)
               cosd(α)  -sind(α)]
    
    u_wing,v_wing = Ṙ * [x_inertial,y_inertial] - R * [U∞,W∞] - Ṙ * [X₀,Y₀]
    return  u_wing,v_wing
end 

function airfoil(naca_digits::String="0012";α=0,chord=1,n_panels=50,x_LE=0,y_LE=0,pitch_axis=0.5,leading_edge_shedding=true)
    n_panels = Int(n_panels)
    n_panels = iseven(n_panels) ? n_panels : n_panels + 1 
 
    x_wing,y_wing = nacacamberline(naca_digits,chord=chord,n_panels=n_panels) # Generate coordinates zero at LE

    x_pitch_axis_wing = chord * pitch_axis # Generate center of rotation in current coordinate system
    y_pitch_axis_wing = 0
 

    x_wing = x_wing .- x_pitch_axis_wing # shift origin to center of rotation
    y_wing = y_wing .- y_pitch_axis_wing

    x_pitch_axis_wing = 0 # Update center of rotation in shifted coordinates
    y_pitch_axis_wing = 0

    x_pitch_axis_inertial, y_pitch_axis_inertial = wing_to_inertial_position(x_pitch_axis_wing,y_pitch_axis_wing,x_LE,y_LE,α) 
    
    
    
    inertial = wing_to_inertial_position.(x_wing,y_wing,x_LE,y_LE,α) 
    x_inertial = [inertial[i][1] for i in 1:length(inertial)]
    y_inertial = [inertial[i][2] for i in 1:length(inertial)]

    all_normals_x_wing = zeros(n_panels)
    all_normals_y_wing = zeros(n_panels)
    all_tangent_x_wing = zeros(n_panels)
    all_tangent_y_wing = zeros(n_panels)
    all_normals_x_inertial = zeros(n_panels)
    all_normals_y_inertial = zeros(n_panels)
    all_tangent_x_inertial = zeros(n_panels)
    all_tangent_y_inertial = zeros(n_panels)


    for i in 1:n_panels-1
        tangentX = x_wing[i+1] - x_wing[i]
        tangentY = y_wing[i+1] - y_wing[i]
        magnitude = sqrt(tangentX^2 + tangentY^2)
        all_tangent_x_wing[i+1] = tangentX/magnitude
        all_tangent_y_wing[i+1] = tangentY/magnitude
        all_normals_x_wing[i+1] = -tangentY/magnitude
        all_normals_y_wing[i+1] = tangentX/magnitude
    end 
    all_tangent_x_wing[1] = all_tangent_x_wing[2]
    all_tangent_y_wing[1] = all_tangent_y_wing[2]
    all_normals_x_wing[1] = all_normals_x_wing[2]
    all_normals_y_wing[1] = all_normals_y_wing[2]


    
    all_normals_inertial = rotate_vector.(all_normals_x_wing,all_normals_y_wing,α)
    all_tangent_inertial = rotate_vector.(all_tangent_x_wing,all_tangent_y_wing,α)

    all_normals_x_inertial = [all_normals_inertial[i][1] for i in 1:length(all_normals_inertial)]
    all_normals_y_inertial = [all_normals_inertial[i][2] for i in 1:length(all_normals_inertial)]

    all_tangent_x_inertial = [all_tangent_inertial[i][1] for i in 1:length(all_tangent_inertial)]
    all_tangent_y_inertial = [all_tangent_inertial[i][2] for i in 1:length(all_tangent_inertial)]
    
    num_vortices = Int((n_panels) / 2) 
    num_cols = Int((n_panels) / 2) 

    empty_vectors_vortex = [Float64[] for i in 1:num_vortices]
    empty_vectors_collocation = [Float64[] for i in 1:num_cols]



    vortices = Vortex.(x_wing[1:2:end-1],
                            y_wing[1:2:end-1],
                            x_inertial[1:2:end-1],
                            y_inertial[1:2:end-1],
                            deepcopy(empty_vectors_vortex),
                            deepcopy(empty_vectors_vortex),
                            all_normals_x_wing[1:2:end-1],
                            all_normals_y_wing[1:2:end-1],
                            all_normals_x_inertial[1:2:end-1],
                            all_normals_y_inertial[1:2:end-1],
                            all_tangent_x_wing[1:2:end-1],
                            all_tangent_y_wing[1:2:end-1],
                            all_tangent_x_inertial[1:2:end-1],
                            all_tangent_y_inertial[1:2:end-1],
                            1.)

    collocations = CollocationPoint.(x_wing[2:2:end],
                            y_wing[2:2:end],
                            x_inertial[2:2:end],
                            y_inertial[2:2:end],
                            deepcopy(empty_vectors_collocation),
                            deepcopy(empty_vectors_collocation),
                            all_normals_x_wing[2:2:end],
                            all_normals_y_wing[2:2:end],
                            all_normals_x_inertial[2:2:end],
                            all_normals_y_inertial[2:2:end],
                            all_tangent_x_wing[2:2:end],
                            all_tangent_y_wing[2:2:end],
                            all_tangent_x_inertial[2:2:end],
                            all_tangent_y_inertial[2:2:end],
                            )
     
    return Airfoil(naca_digits,chord,pitch_axis,x_pitch_axis_inertial,y_pitch_axis_inertial,vortices,collocations,leading_edge_shedding)
end  

function rotate_vector(x,y,α)
    R = [cosd(α) sind(α)
        -sind(α) cosd(α)]

    x_rot,y_rot = R * [x,y]
    return x_rot,y_rot
end 

function draw(airfoil::Airfoil;inertial_frame=true)
    num = inertial_frame ? 2 : 0
    x_vort = [getfield(vort,1+num) for vort in airfoil.vortices]
    y_vort = [getfield(vort,2+num) for vort in airfoil.vortices]
    x_col = [getfield(col,1+num) for col in airfoil.collocations]
    y_col = [getfield(col,2+num) for col in airfoil.collocations]
    h = scatter(x_vort,y_vort,aspect_ratio=:equal,label="Bound Vortices")
    scatter!(x_col,y_col,aspect_ratio=:equal,label="Collocation Points")
    display(h)
end 

function draw(flow::Flowfield;inertial_frame=true)
    num = inertial_frame ? 2 : 0
    x_vort = [getfield(vort,1+num) for vort in flow.wing.vortices]
    y_vort = [getfield(vort,2+num) for vort in flow.wing.vortices]
    x_col = [getfield(col,1+num) for col in flow.wing.collocations]
    y_col = [getfield(col,2+num) for col in flow.wing.collocations]
    x_wake_vortices = [getfield(vort,1+num) for vort in flow.wake_vortices]
    y_wake_vortices = [getfield(vort,2+num) for vort in flow.wake_vortices]
    h = scatter(x_vort,y_vort,aspect_ratio=:equal,label="Bound Vortices")
    scatter!(x_col,y_col,aspect_ratio=:equal,label="Collocation Points")
    scatter!(x_wake_vortices,y_wake_vortices,aspect_ratio=:equal,label="Shed Vortices")
    display(h)
end 

function induced_velocity(target,source::Vortex;unit_forcing=true,inertial_frame=false)
    A = [0 1
        -1 0]
    
    num = inertial_frame ? 2 : 0
    target_x = getfield(target,1+num)
    target_y = getfield(target,2+num)
    source_x = getfield(source,1+num)
    source_y = getfield(source,2+num)

    Γ = unit_forcing ? 1 : source.Γ
    X = [target_x - source_x ,target_y - source_y]
    r² = (target_x - source_x)^2 + (target_y - source_y)^2
    C = Γ/(2π*r²)
    u,v = C*A*X
    return u,v
end 

function movewing!(flow::Flowfield,kinematics::Kinematics,nstep::Int64)
    α = kinematics.α[nstep]
    
    flow.wing.x_pitch_axis = kinematics.X₀[nstep]
    flow.wing.y_pitch_axis = kinematics.Y₀[nstep]

    for vort in flow.wing.vortices
        vort.x_inertial,vort.y_inertial = wing_to_inertial_position(vort.x_wing,vort.y_wing,flow.wing.x_pitch_axis,flow.wing.y_pitch_axis,α)
        vort.x_normal_inertial,vort.y_normal_inertial = rotate_vector(vort.x_normal_wing,vort.y_normal_wing,α)
        vort.x_tangent_inertial,vort.y_tangent_inertial = rotate_vector(vort.x_tangent_wing,vort.y_tangent_wing,α)
        push!(vort.x_inertial_history,vort.x_inertial)
        push!(vort.y_inertial_history,vort.y_inertial)
    
    end 

    for col in flow.wing.collocations
        col.x_inertial,col.y_inertial = wing_to_inertial_position(col.x_wing,col.y_wing,flow.wing.x_pitch_axis,flow.wing.y_pitch_axis,α)
        col.x_normal_inertial,col.y_normal_inertial = rotate_vector(col.x_normal_wing,col.y_normal_wing,α)
        col.x_tangent_inertial,col.y_tangent_inertial = rotate_vector(col.x_tangent_wing,col.y_tangent_wing,α)
        push!(col.x_inertial_history,col.x_inertial)
        push!(col.y_inertial_history,col.y_inertial)
    end 
end 


function placevortex!(flow::Flowfield,kinematics::Kinematics,nstep)
    TE_x_travel = flow.wing.collocations[end].x_inertial_history[end] - flow.wing.collocations[end].x_inertial_history[end-1]
    TE_y_travel = flow.wing.collocations[end].y_inertial_history[end] - flow.wing.collocations[end].y_inertial_history[end-1]
    ratio = 0.7 
    new_vort_x_inertial = ratio*TE_x_travel + flow.wing.collocations[end].x_inertial_history[end-1]
    new_vort_y_inertial = ratio*TE_y_travel + flow.wing.collocations[end].y_inertial_history[end-1]
    new_vort_x_wing,new_vort_y_wing = inertial_to_wing_position(new_vort_x_inertial,new_vort_y_inertial,kinematics.X₀[nstep],kinematics.Y₀[nstep],kinematics.α[nstep])
    push!(flow.wake_vortices,Vortex(new_vort_x_wing,
                             new_vort_y_wing,
                             new_vort_x_inertial,
                             new_vort_y_inertial,
                             Float64[new_vort_x_inertial],
                             Float64[new_vort_y_inertial],
                             0.,
                             0.,
                             0.,
                             0.,
                             0.,
                             0.,
                             0.,
                             0.,
                             1.)
    )
end 


function enforce_no_throughflow!(flow::Flowfield,kinematics::Kinematics,nstep)
    A = populate_influence_matrix(flow)
    RHS = populate_RHS(flow,kinematics,nstep)
    solve_system!(A,RHS,flow)
end 

function populate_influence_matrix(flow::Flowfield)
    A = zeros(length(flow.wing.vortices)+1,length(flow.wing.vortices)+1)
    for (i,col) in enumerate(flow.wing.collocations)
        for (j,vort) in enumerate(flow.wing.vortices)
            u,v = induced_velocity(col,vort)
            A[i,j] = dot([u,v],[col.x_normal_wing,col.y_normal_wing])
        end 
    end 

    for (i,col) in enumerate(flow.wing.collocations)
        u,v = induced_velocity(col,flow.wake_vortices[end])
        A[i,end] = dot([u,v],[col.x_normal_wing,col.y_normal_wing])
    end 
    A[end,:] .= 1
    return A
end 

function populate_RHS(flow::Flowfield,kinematics::Kinematics,nstep)
    RHS = zeros(length(flow.wing.vortices)+1)
    α = kinematics.α[nstep]
    R = [cosd(α) -sind(α)
         sind(α)  cosd(α)]

    for (i,col) in enumerate(flow.wing.collocations)
        uₖ,vₖ = R*[-kinematics.U₀[nstep],-kinematics.V₀[nstep]] + [-kinematics.α̇[nstep]*col.y_wing,kinematics.α̇[nstep]*col.x_wing]
        uᵥ = 0
        vᵥ = 0
        for vort in flow.wake_vortices[1:end-1]
            u_vort,v_vort = induced_velocity(col,vort,unit_forcing=false)
            uᵥ += u_vort
            vᵥ += v_vort
        end 
        RHS[i] = -dot([uᵥ+uₖ,vᵥ+vₖ],[col.x_normal_wing,col.y_normal_wing])
    end 

    RHS[end] = (length(flow.wake_vortices)-1) ==0 ? 0 : -sum([flow.wake_vortices[k].Γ for k in 1:length(flow.wake_vortices)-1])

    return RHS
end

function solve_system!(A::Array{Float64,2},RHS::Array{Float64,1},flow::Flowfield)
    Γs = A \ RHS
    for (i,vort) in enumerate(flow.wing.vortices)
        vort.Γ = Γs[i]
    end 
    flow.wake_vortices[end].Γ = Γs[end]
end 

function wake_rollup!(flow::Flowfield,kinematcs::Kinematics,Δt::Float64,nstep)
    all_u = zeros(length(flow.wake_vortices))
    all_v = zeros(length(flow.wake_vortices))
    for (i,wakevort) in enumerate(flow.wake_vortices)
        u_total = 0
        v_total = 0
        for (j,vort) in enumerate(flow.wake_vortices)    
            if i != j
                u_induced, v_induced= induced_velocity(wakevort,vort,inertial_frame=true,unit_forcing=false)
                u_total += u_induced
                v_total += v_induced
            end 
        end

        for vort in flow.wing.vortices
            u_induced, v_induced= induced_velocity(wakevort,vort,inertial_frame=true,unit_forcing=false)
            u_total += u_induced
            v_total += v_induced
        end 

        all_u[i] = u_total
        all_v[i] = v_total
    end 

    for (i,wakevort) in enumerate(flow.wake_vortices)
        wakevort.x_inertial += Δt * all_u[i]
        wakevort.y_inertial += Δt * all_v[i]
        append!(wakevort.x_inertial_history,wakevort.x_inertial)
        append!(wakevort.y_inertial_history,wakevort.y_inertial)

        x_wing,y_wing = inertial_to_wing_position(wakevort.x_inertial, wakevort.y_inertial,flow.wing.x_pitch_axis,flow.wing.y_pitch_axis,kinematcs.α[nstep])
        wakevort.x_wing = x_wing
        wakevort.y_wing = y_wing
    end 


    if flow.wing.leading_edge_shedding
        LEV = deepcopy(flow.wing.vortices[1])
        u_total = 0
        v_total = 0
        for (i,wakevort) in enumerate(flow.wake_vortices)
            u_induced, v_induced= induced_velocity(LEV,wakevort,inertial_frame=true,unit_forcing=false)
            u_total += u_induced
            v_total += v_induced 
        end
    
        for vort in flow.wing.vortices[2:end]
            u_induced, v_induced= induced_velocity(LEV,vort,inertial_frame=true,unit_forcing=false)
            u_total += u_induced
            v_total += v_induced
        end 

        LEV.x_inertial += Δt * u_total
        LEV.y_inertial += Δt * v_total
        append!(LEV.x_inertial_history,LEV.x_inertial)
        append!(LEV.y_inertial_history,LEV.y_inertial)

        x_wing,y_wing = inertial_to_wing_position(LEV.x_inertial, LEV.y_inertial,flow.wing.x_pitch_axis,flow.wing.y_pitch_axis,kinematcs.α[nstep])
        LEV.x_wing = x_wing
        LEV.y_wing = y_wing

        push!(flow.wake_vortices,LEV)

    end 


end

function runme()
    foil = airfoil(leading_edge_shedding=true,α=30) 
    Δt = 0.01
    time_steps  = 300
    final_time = time_steps * Δt
    mytime = collect(0:Δt:15)
    AoA = 20*sin.(10*mytime)#30*ones(length(mytime))#
    kinematics = set_kinematics(mytime,-1,1,AoA,foil)
    gust = top_hat(kinematics.X₀)
    flow = Flowfield(foil,Vortex[],gust)
    movewing!(flow,kinematics,1)
    draw(flow)
    for nstep in 2:time_steps
        movewing!(flow,kinematics,nstep)
        placevortex!(flow,kinematics,nstep)
        draw(flow,inertial_frame=true)
        #enforce_no_throughflow!(flow,kinematics,nstep)
        wake_rollup!(flow,kinematics,Δt,nstep)
       # draw(flow,inertial_frame=true)
    end 
    return flow,kinematics
end  


# Possible bugs