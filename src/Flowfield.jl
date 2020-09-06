module Flowfield

include("./Util.jl")
using .Util

using Makie, AbstractPlotting

using OffsetArrays, CoordinateTransformations, LinearAlgebra
using Random
using Plots
using Luxor
using Images, Colors, ImageFiltering
using Interpolations


function generateRain(heightmap, seed)
        rainThreshold = 0.7
        HeightBias = 10
        mapSize = size(heightmap)

        scale = maximum(heightmap) - minimum(heightmap)
        scaled = heightmap / (HeightBias*scale)

        return clamp.(sqrt.(clamp.(
        Util.generateNoise(mapSize[1], mapSize[2], seed; freq=0.1, d=4, pers=0.5) .- rainThreshold,
        0, 1)), 0,1)
end

function erosionDepthRamp(x, Kdmax)

        return clamp.(1 .- x/Kdmax, 0,1)


        if x <= 0
                return 0
        elseif x >= Kdmax
                return 1
        else
                1 .- (Kdmax - x)./Kdmax
        end

end



function erodeTerrain(map, iterations=100, resetSeed=false)
        #Initialize (mapSize, resetSeed)
        # http://www-ljk.imag.fr/Publications/Basilic/com.lmc.publi.PUBLI_Inproceedings@117681e94b6_fff75c/FastErosion_PG07.pdf

        mapSize = size(map)

        # parameters
        dt::Float32 = 0.005
        g::Float32 = 9.81
        crossSection::Float32 = 5
        sedimentCapacity::Float32 = 0.5
        minAngle::Float32 = 0.2
        sinMinAngle = sin(minAngle)
        erosionFactor::Float32 = 0.5
        depositionFactor::Float32 = 1
        evaporationFactor::Float32 = 0.2
        rainThreshold::Float32 = 0.7
        rainIntensity::Float32 = 0.1
        thermalErosion::Float32 = 0.1
        inertia::Float32 = 0.5
        maxSoftness::Float32 = 5
        softening::Float32 = 0.5

        KdMax::Float32 = 20
        InvKdMax = 1/KdMax


        talusAngle::Float32 = 70*pi/180
        talusHeight::Float32 = tan(talusAngle)
        talusAngleBias::Float32 = 0.1

        neigbourCells = [(-1, -1), (0, -1), (1, -1), (-1, 1), (0, 1), (1, 1), (1, 0), (-1, 0)]
        neigbourCellsInverse = [(1, 1), (0, 1), (-1, 1), (1, -1), (0, -1), (-1, -1), (-1, 0), (1, 0)]

        # calculation constants

        maxVel::Float32 = 1/dt

        height = copy(map)
        waterHeight = zeros(Float32, mapSize)
        sediment = zeros(Float32, mapSize)
        maxCapacity = zeros(Float32, mapSize)

        f_L = zeros(Float32, mapSize)
        f_R = zeros(Float32, mapSize)
        f_T = zeros(Float32, mapSize)
        f_B = zeros(Float32, mapSize)
        dh_L = zeros(Float32, mapSize)
        dh_R = zeros(Float32, mapSize)
        dh_T = zeros(Float32, mapSize)
        dh_B = zeros(Float32, mapSize)
        flowScaleK = zeros(Float32, mapSize)

        velX = zeros(Float32, mapSize)
        velY = zeros(Float32, mapSize)
        velEro = zeros(Float32, mapSize)

        softness = ones(Float32, mapSize)

        #heightPlot = Node(height)
        #waterHeightPlot = Node(waterHeight)
        scene = Makie.Scene(resolution=(1500, 1200))
        waterInflux = rainIntensity*generateRain(height, 123123)


        for iter = 1:iterations
                #println("iter ", iter, "/", iterations)

                #println("-------------------------")


                if iter % 10 == 0
                        println("iter ", iter, "/", iterations)
                end


                if 50 < iter < iterations - 10
                        waterInflux = rainIntensity*ones(mapSize)
                #elseif iter == iterations - 100
                #        waterInflux = 100*rainIntensity*ones(mapSize)
                elseif iter > iterations - 10
                        waterInflux = zeros(mapSize)
                        #evaporationFactor = 0
                        #erosionFactor = 0
                end

                #elseif iter % 5 == 1
                #        waterInflux = rainIntensity*generateRain(height, iter)
                #else
                #        waterInflux = zeros(mapSize[1], mapSize[2])
                #end
                # add water

                waterHeight1 = waterHeight + waterInflux*dt

                #if iter % 10 == 0
                #        plotMesh(height, waterHeight1, sediment)
                #end
                #display(waterHeight1)

                #println("max(waterHeight1) = ", maximum(abs.(waterHeight1)))

                # calculate slopes
                dh_L[:, 2:end] =   height[:, 2:end]   + waterHeight1[:, 2:end]   - height[:, 1:end-1] - waterHeight1[:, 1:end-1]
                dh_R[:, 1:end-1] = height[:, 1:end-1] + waterHeight1[:, 1:end-1] - height[:, 2:end]   - waterHeight1[:, 2:end]
                dh_T[2:end, :] =   height[2:end, :]   + waterHeight1[2:end, :]   - height[1:end-1, :] - waterHeight1[1:end-1, :]
                dh_B[1:end-1, :] = height[1:end-1, :] + waterHeight1[1:end-1, :] - height[2:end, :]   - waterHeight1[2:end, :]


                fact = dt*crossSection*g

                # calculate directional flow rates
                f_L = max.(0, f_L + fact*dh_L)
                f_R = max.(0, f_R + fact*dh_R)
                f_T = max.(0, f_T + fact*dh_T)
                f_B = max.(0, f_B + fact*dh_B)

                #println("max(f_L) = ", maximum(abs.(f_L)))

                # scale flow rate depending on total flow
                totalFlow = dt*(f_L + f_R + f_T + f_B)

                #println(size(totalFlow))
                #println(totalFlow)
                #println("min(totalFlow) = ", minimum(abs.(totalFlow)))

                for idx in eachindex(flowScaleK)
                        if totalFlow[idx] < waterHeight1[idx] || totalFlow[idx] == 0
                                flowScaleK[idx] = 1
                        else
                                flowScaleK[idx] = waterHeight1[idx]/totalFlow[idx]
                        end
                end


                f_L = f_L .* flowScaleK
                f_R = f_R .* flowScaleK
                f_T = f_T .* flowScaleK
                f_B = f_B .* flowScaleK

                #println(maximum(abs.(waterHeight1)))

                # shifted flux for summing
                f_L_roll = circshift(f_L, (0, -1))
                f_R_roll = circshift(f_R, (0, 1))
                f_T_roll = circshift(f_T, (-1, 0))
                f_B_roll = circshift(f_B, (1, 0))

                # calculate volumetric flow sum of inward flow from adjacent cells minus outward flow



                inflow = (f_L_roll + f_R_roll + f_T_roll + f_B_roll)
                outflow = (f_L + f_R + f_T + f_B)

                dV = dt*(inflow - outflow)

                #dV = inflow - outflow
                waterHeight2 = max.(0, waterHeight1 + dV)
                #display(waterHeight2)

                #display(waterHeight2)


                #println(waterHeight2)

                # update valocity
                avgWaterHeight = 0.5*(waterHeight2 + waterHeight1)

                dWX = 0.5*(f_R_roll + f_R - f_L - f_L_roll)
                dWY = 0.5*(f_B_roll + f_B - f_T - f_T_roll)

                #println("max(dwx) = ", maximum(abs.(dWX)))
                #println("max(dwx) = ", maximum(abs.(dWX)))
                #println("max(avgWaterHeight) = ", maximum(abs.(avgWaterHeight)))

                for idx in eachindex(velX)
                        if dWX[idx] == 0 || avgWaterHeight[idx] == 0
                                velX[idx] = 0
                        else
                                velX[idx] =  clamp(dWX[idx] / avgWaterHeight[idx], -maxVel, maxVel)
                        end
                end

                for idx in eachindex(velY)
                        if dWY[idx] == 0 || avgWaterHeight[idx] == 0
                                velY[idx] = 0
                        else
                                velY[idx] =  clamp(dWY[idx] / avgWaterHeight[idx], -maxVel, maxVel)
                        end
                end

                for idx in eachindex(velEro)
                        if outflow[idx] == 0 || avgWaterHeight[idx] == 0
                                velEro[idx] = 0
                        else
                                velEro[idx] =  clamp(outflow[idx] / avgWaterHeight[idx], -maxVel, maxVel)
                        end
                end


                #lmax = max.(1 .- waterHeight1*InvKdMax,1)
                # calculate sediment capacitymax.(1 .- waterHeight1*InvKdMax,1)
                #maxCapacity = sedimentCapacity .* max.(Util.sinLocalSlope(height), sinMinAngle) .* sqrt.(velX.^2 + velY.^2)
                maxCapacity = sedimentCapacity .* max.(Util.sinLocalSlope(height), sinMinAngle) .*
                        max.(velEro, sqrt.(velX.^2 + velY.^2))# .* lmax


                vel = zeros(Float32, mapSize[1], mapSize[2], 3)
                vel[:,:,1] = velX
                vel[:,:,2] = velY

                neg_norm = Util.normal(height)
                vel_3d = Util.projectionOnPlane.(vel, neg_norm)
                coll_3d = dot.(-Util.normal.(height), Util.projectionOnPlane.(vel))

                maxCapacity = sedimentCapacity .* coll_3d .*
                        max.(velEro, sqrt.(velX.^2 + velY.^2))# .* lmax


                capacityDiff = (sediment-maxCapacity)
                # deposit and erode
                depositAmount = ifelse.(capacityDiff .> 0,
                        dt*depositionFactor.*softness.*capacityDiff,
                        dt*erosionFactor*capacityDiff
                )

                depositAmount = min.(depositAmount, waterHeight2)


                height += depositAmount
                sediment1 = sediment - depositAmount
                waterHeight3 = waterHeight2 - depositAmount

                # increase softness when sediment is deposited
                softness = ifelse.(capacityDiff .> 0,
                        max.(maxSoftness, softness  + dt*erosionFactor.*softening.*capacityDiff),
                        softness
                )


                talusFlowBuffer = zeros(Float32, mapSize[1], mapSize[2], 8)

                # calculate sediment transport
                sediment1Interp = interpolate(sediment1, BSpline(Linear()))
                for (i,j) in Iterators.product(1:mapSize[1], 1:mapSize[2])
                        sediment[i,j] = sediment1Interp(
                        clamp(i - (velX[i,j] * dt), 1, size(sediment, 1)),
                        clamp(j - (velY[i,j] * dt), 1, size(sediment, 2))
                        )

                        # thermal erosion by talus angle

                        # max_dh_talus = 0
                        # for k in 1:9
                        #         ik = i + neigbourCells[k][1]
                        #         jk = j + neigbourCells[k][2]
                        #         if 1 <= ik <= mapSize[1] && 1 <= jk <= mapSize[2]
                        #                 dh_talus = height[i,j] - height[ik, jk]
                        #                 if dh_talus > norm(neigbourCells[k]) * talusHeight
                        #                         talusFlowBuffer[i,j,k] = dh_talus
                        #                 end
                        #                 if dh_talus > max_dh_talus
                        #                         max_dh_talus = dh_talus
                        #                 end
                        #         end
                        # end
                        #
                        # talusVolumeFlow = 0.5*dt*thermalErosion*maxiumum(talusFlowBuffer[i,j,:])
                        # talusFlowSum = sum(talusFlowBuffer[i,j,:])
                        # if talusFlowSum > 0
                        #         talusFlowBuffer[i,j,:] = talusVolumeFlow .* talusFlowBuffer[i,j,:] ./ talusFlowSum
                        # end
                end




                # evaporate water
                waterHeight = waterHeight3 * (1 - (evaporationFactor*dt))


                #display(waterHeight)
                #println(sum(waterHeight))

        end

        height += sediment
        waterHeight -= sediment

        println(minimum(waterHeight))
        println(maximum(height))
        plotMesh(height, waterHeight, sediment)
        #plotHeightAndWaterInit(height, waterHeight, mapSize[1], mapSize[2], scene)

        return height, waterHeight
end

function plotHeightAndWaterUpdate(heightmap, water, sx, sy, scene)

end

function plotHeightAndWaterInit(heightmap, water, sx, sy, scene)

        ambient = 0.4
        diffuse = 0.6
        specular = 0.05
        shininess = 8.0

        radius = 3.0
        theta = 30.0
        phi = 90.0

        la = map(Vec3f0, ambient)
        ld = map(Vec3f0, diffuse)
        ls = map(Vec3f0, specular)
        lp = map(radius, theta, phi) do r, theta, phi
                r * Vec3f0(
                cosd(phi) * sind(theta),
                sind(phi) * sind(theta),
                cosd(theta)
                )
        end

        scene2 = AbstractPlotting.surface(
        range(0,1,length=sx), range(0,1,length=sy), (heightmap)./sx,
        ambient = la, diffuse = ld, specular = ls, shininess = shininess,
        lightposition = lp,
        shading=true,
        #color=fill(RGB(0.4,0.4,1.0),sx,sy)
        #showaxes=false,
        #show_axis=false
        )

        # AbstractPlotting.surface!(scene2,
        #     range(0,1,length=sx), range(0,1,length=sy), heightmap./sx,
        #     ambient = la, diffuse = ld, specular = ls, shininess = shininess,
        #     lightposition = lp,
        #     shading=true,
        #     color = fill(RGBA(1.,1.,1.,1),sx,sy),
        #     transparency=false
        #     #showaxes=false,
        #     #show_axis=false
        # )


        #Makie.scatter!(scene2, map(v -> [v], lp), color=:yellow, markersize=1f0)

        # Combine scene
        vbox(hbox(scene2), parent=scene)

        AbstractPlotting.inline!(false)
        display(scene)
        println("End plot")


        #scene
        # Do not execute beyond this point!
        #RecordEvents(scene, "lighting")
        println("done plotting")
end

function plotHeightmap(heightmap, sx, sy, transparent=false)
        println("start plotting")
        # Set up sliders to control lighting attributes
        s1, ambient = textslider(0f0:0.01f0:1f0, "ambient", start = 0.4f0)
        s2, diffuse = textslider(0f0:0.025f0:2f0, "diffuse", start = 0.6f0)
        s3, specular = textslider(0f0:0.025f0:2f0, "specular", start = 0.05f0)
        s4, shininess = textslider(2f0.^(2f0:8f0), "shininess", start = 8f0)

        # Set up (r, θ, ϕ) for lightposition
        s5, radius = textslider(2f0.^(0.5f0:0.25f0:20f0), "light pos r", start = 3f0)
        s6, theta = textslider(0:5:180, "light pos theta", start = 60f0)
        s7, phi = textslider(0:5:360, "light pos phi", start = 70f0)

        # transform signals into required types
        la = map(Vec3f0, ambient)
        ld = map(Vec3f0, diffuse)
        ls = map(Vec3f0, specular)
        lp = map(radius, theta, phi) do r, theta, phi
                r * Vec3f0(
                cosd(phi) * sind(theta),
                sind(phi) * sind(theta),
                cosd(theta)
                )
        end

        println("done plotting1")

        # Set up surface plot + light source
        scene2 = AbstractPlotting.surface(
        range(0,1,length=sx), range(0,1,length=sy), heightmap,
        ambient = la, diffuse = ld, specular = ls, shininess = shininess,
        lightposition = lp,
        shading=true,
        color=Gray,
        transparency=transparent,
        alpha=0.5
        #showaxes=false,
        #show_axis=false
        )
        #Makie.scatter!(scene2, map(v -> [v], lp), color=:yellow, markersize=1f0)

        # Combine scene
        scene = Makie.Scene(resolution=(1500, 1200))
        vbox(hbox(s4, s3, s2, s1, s7, s6, s5), hbox(scene2), parent=scene)

        AbstractPlotting.inline!(false)
        display(scene)
        println("End plot")


        #scene
        # Do not execute beyond this point!
        RecordEvents(scene, "lighting")
        println("done plotting")
end

function plotHeightmap2(heightmap, sx, sy)
        println("Begin plot")


        p1 = AbstractPlotting.surface(
        1:sx, 1:sy, heightmap,
        colormap = :Spectral
        )
        scene = vbox(p1)
        #popdisplay(AbstractPlotting.PlotDisplay())
        AbstractPlotting.inline!(false)
        display(scene)
        println("End plot")


        scene

end

function plotMesh(heightmap, waterHeight, sediment)

        meshedPts, connectPts = Util.generateMesh(heightmap, 1)

        # Set up sliders to control lighting attributes
        s1, ambient = textslider(0f0:0.01f0:1f0, "ambient", start = 0.55f0)
        s2, diffuse = textslider(0f0:0.025f0:2f0, "diffuse", start = 0.6f0)
        s3, specular = textslider(0f0:0.025f0:2f0, "specular", start = 0.2f0)
        s4, shininess = textslider(2f0.^(2f0:8f0), "shininess", start = 8f0)

        # Set up (r, θ, ϕ) for lightposition
        s5, radius = textslider(2f0.^(0.5f0:0.25f0:20f0), "light pos r", start = 5000f0)
        s6, theta = textslider(0:5:180, "light pos theta", start = 45f0)
        s7, phi = textslider(0:5:360, "light pos phi", start = 65f0)

        # transform signals into required types
        la = map(Vec3f0, ambient)
        ld = map(Vec3f0, diffuse)
        ls = map(Vec3f0, specular)
        lp = map(radius, theta, phi) do r, theta, phi
                r * Vec3f0(
                cosd(phi) * sind(theta),
                sind(phi) * sind(theta),
                cosd(theta)
                )
        end

        scene2 = AbstractPlotting.mesh(meshedPts, connectPts,
        shading = true,
        color=:gray,
        ambient = la, diffuse = ld, specular = ls, shininess = shininess,
        lightposition = lp,
        showaxes=false,
        show_axis=false
        )



        scale = maximum(heightmap) - minimum(heightmap)

        if true
                waterLevel = minimum(heightmap).+15
                waterPlot = waterLevel.*ones(size(heightmap))

                waterColorPlot = max.(waterLevel .- heightmap, 0)

                meshedPtsWOnly, _ = Util.generateMesh(max.(waterHeight,0), 1)

                color = meshedPtsWOnly[:, 3] ./ maximum(meshedPtsWOnly[:, 3])

        else
                waterPlot = ifelse.(waterHeight .< eps(Float32),
                0, heightmap + waterHeight)

                meshedPtsWOnly, _ = Util.generateMesh(max.(waterHeight,0), 1)

                color = meshedPtsWOnly[:, 3] ./ maximum(meshedPtsWOnly[:, 3])
        end

        meshedPtsW, connectPtsW = Util.generateMesh(waterPlot, 1)



        AbstractPlotting.mesh!(scene2,
        meshedPtsW, connectPtsW,
        shading = true,
        color=color,
        colormap=[RGBA(0.18, 0.45, 0.62, 0.5), RGBA(0.18, 0.45, 0.62, 0.95)],
        ambient = la, diffuse = ld, specular = ls, shininess = shininess,
        lightposition = lp,
        showaxes=false,
        show_axis=false
        )

        # Combine scene
        scene = Makie.Scene(resolution=(1500, 1200))
        vbox(hbox(scene2), parent=scene)

        AbstractPlotting.inline!(false)
        display(scene)
        println("End plot")


        #scene
        # Do not execute beyond this point!
        RecordEvents(scene, "lighting")
        println("done plotting")

end

println("Begin ====================================================")

sx = 256
sy = 256
hScale = sx/2
seed = 123135
numIterations = 400
Random.seed!(seed)

#rainThreshold = 0.75
#waterInflux = sqrt.(clamp.(Util.generateNoise(sx, sy, seed+6; freq=0.05, d=4, pers=0.99) .- rainThreshold, 0, 1))
#plotHeightmap(waterInflux, sx, sy)


#erosionBrush = initializeBrush(erosionRadius)
#depositionBrush = initializeBrush(erosionRadius)
heightmap = hScale*Util.generateNoise(sx, sy, seed, freq=0.008, d=6, pers=0.45)

#plotMesh(heightmap, zeros(size(heightmap)), zeros(size(heightmap)))
heightmapEroded, waterHeight = erodeTerrain(heightmap, numIterations)

change = heightmapEroded - heightmap

#plotMesh(change, zeros(sx,sy), zeros(sx,sy))


#diffHeight = image2 - heightmap
#plot(image)


#plotHeightmap(heightmapEroded./sx, sx, sy, false)
#plotHeightmap(2*image2./sx, sx, sy, false)
#plotMesh(image2)
println("done")


end
