module ParticleSinglestep

using Makie, AbstractPlotting

using OffsetArrays, CoordinateTransformations
using Random
using Plots
using Luxor
using Images, Colors, ImageFiltering

struct HeightAndGradient
        height::AbstractFloat
        gradientX::AbstractFloat
        gradientY::AbstractFloat
end


erosionRadius = 6
depositionRadius = 4
inertia = 0.2 # At zero, water will instantly change direction to flow downhill. At 1, water will never change direction.
sedimentCapacityFactor = 32.0 # Multiplier for how much sediment a droplet can carry
minSlope = 0.01 # Used to prevent carry capacity getting too close to zero on flatter terrain

erodeSpeed = 1
depositSpeed = 0.6
evaporateSpeed = 0.04
gravity = 10.0
maxDropletLifetime = 256
friction = 0.0

initialWaterVolume = 20.0
initialSpeed = 3
speedThreshold = 0.001
blurFactor = 0.2
blurRadius = 5.0






function calcHeightAndGradient(heightmap, posX::AbstractFloat, posY::AbstractFloat)::HeightAndGradient

        float_type = typeof(heightmap[1,1])

        coordX::Int = floor(posX)
        coordY::Int = floor(posY)

        # Calculate droplet's offset inside the cell (0,0) = at NW node, (1,1) = at SE node
        x = posX - coordX
        y = posY - coordY

        # Calculate heights of the four nodes of the droplet's cell
        heightNW = heightmap[coordX, coordY]
        heightNE = heightmap[coordX+1, coordY]
        heightSW = heightmap[coordX, coordY+1]
        heightSE = heightmap[coordX+1, coordY+1]

        # Calculate droplet's direction of flow with bilinear interpolation of height difference along the edges
        gradientX::float_type = (heightNE - heightNW) * (1 - y) + (heightSE - heightSW) * y;
        gradientY::float_type = (heightSW - heightNW) * (1 - x) + (heightSE - heightNE) * x;

        # Calculate height with bilinear interpolation of the heights of the nodes of the cell
        height::float_type = heightNW * (1 - x) * (1 - y) + heightNE * x * (1 - y) + heightSW * (1 - x) * y + heightSE * x * y;

        return HeightAndGradient(height, gradientX, gradientY)

end

function initializeBrush(radius)

    radiusRange = -radius+1:radius-1

    function circleWeight(dx, dy)
        return max(0, 1 - sqrt(dx^2 + dy^2) / radius)
    end

    wCalc = [circleWeight(x,y) for x in radiusRange, y in radiusRange]
    weightSum = sum(wCalc)
    wCalc /= weightSum
    return OffsetArray(wCalc, radiusRange, radiusRange)

end

function erodeTerrain(map, numIterations=1, resetSeed=false)
        #Initialize (mapSize, resetSeed)
        mapEroded = copy(map)
        mapSize = size(mapEroded)

        for i = 1:numIterations

            if i%1000 == 0
                println("iter ",i, "/", numIterations)
            end


            # Create water droplet at random point on map
            posX = (mapSize[1] - 1)*rand() + 1
            posY = (mapSize[2] - 1)*rand() + 1
            dirX = 0.0
            dirY = 0.0
            speed = initialSpeed
            water = initialWaterVolume
            sediment = 0.0

            #path = Matrix{Float64}(undef, 0, 3)

            for lifetime = 1:maxDropletLifetime
                # if speed < speedThreshold && lifetime < maxDropletLifetime
                #     continue
                # end
                nodeX::Int = floor(posX)
                nodeY::Int = floor(posY)
                dropletIndex = [nodeX nodeY]

                # Calculate droplet's offset inside the cell (0,0) = at NW node, (1,1) = at SE node
                cellOffsetX = posX - nodeX
                cellOffsetY = posY - nodeY

                # Calculate droplet's height and direction of flow with bilinear interpolation of surrounding heights
                heightAndGradient::HeightAndGradient = calcHeightAndGradient(mapEroded, posX, posY)

                #p = [posX  posY  heightAndGradient.height]
                #path = [path;p]

                # Update the droplet's direction and position (move position 1 unit regardless of speed)
                dirX = (dirX * inertia - heightAndGradient.gradientX * (1 - inertia))
                dirY = (dirY * inertia - heightAndGradient.gradientY * (1 - inertia))
                # Normalize direction
                len = sqrt(dirX * dirX + dirY * dirY)
                if (len != 0)
                    dirX /= len
                    dirY /= len
                end
                posX += dirX
                posY += dirY

                # Stop simulating droplet if it's not moving or has flowed over edge of map
                if ((dirX == 0 && dirY == 0) || posX < 1 || posX >= mapSize[1] || posY < 1 || posY >= mapSize[2])
                    break
                end

                # Find the droplet's new height and calculate the deltaHeight
                newHeight = calcHeightAndGradient(mapEroded, posX, posY).height
                deltaHeight = newHeight - heightAndGradient.height;
                #println(heightAndGradient.height)

                # Calculate the droplet's sediment capacity (higher when moving fast down a slope and contains lots of water)
                sedimentCapacity = max(-deltaHeight, minSlope) * speed * water * sedimentCapacityFactor

                # If carrying more sediment than capacity, or if flowing uphill:
                if (sediment > sedimentCapacity || deltaHeight > 0 )
                    # If moving uphill (deltaHeight > 0) try fill up to the current height, otherwise deposit a fraction of the excess sediment
                    if deltaHeight > 0
                        amountToDeposit = min(deltaHeight, sediment)
                        sediment -= amountToDeposit

                        # Add the sediment to the four nodes of the current cell using bilinear interpolation
                        # Deposition is not distributed over a radius (like erosion) so that it can fill small pits
                        mapEroded[dropletIndex[1], dropletIndex[2]] += amountToDeposit * (1 - cellOffsetX) * (1 - cellOffsetY)
                        mapEroded[dropletIndex[1] + 1, dropletIndex[2]] += amountToDeposit * cellOffsetX * (1 - cellOffsetY)
                        mapEroded[dropletIndex[1], dropletIndex[2] + 1] += amountToDeposit * (1 - cellOffsetX) * cellOffsetY
                        mapEroded[dropletIndex[1] + 1, dropletIndex[2] + 1] += amountToDeposit * cellOffsetX * cellOffsetY
                    elseif erosionRadius <= dropletIndex[1] <= sx - erosionRadius && erosionRadius <= dropletIndex[2] <= sy - erosionRadius
                        amountToDeposit = (sediment - sedimentCapacity) * depositSpeed
                        sediment -= amountToDeposit

                        brushTranslated = warp(copy(depositionBrush), Translation(-dropletIndex[1],-dropletIndex[2]))
                        #summary(brushTranslated)
                        brushX = axes(brushTranslated, 1)
                        brushY = axes(brushTranslated, 2)
                        #println(brushX)
                        affectedMap = mapEroded[brushX, brushY]
                        #println(affectedMap)
                        afterDepostionAffected = affectedMap + (amountToDeposit * brushTranslated)

                        # clamp erosion to not go below 0
                        afterDepostionAffected = max.(0, afterDepostionAffected)

                        sediment += sum(affectedMap-afterDepostionAffected)
                        mapEroded[brushX, brushY] = afterDepostionAffected


                    end

                    #mountToDeposit = 0

                    # If last iteration deposit everything
                    # if lifetime == maxDropletLifetime
                    #     amountToDeposit = sediment
                    # end



                else
                    # Erode a fraction of the droplet's current carry capacity.
                    # Clamp the erosion to the change in height so that it doesn't dig a hole in the terrain behind the droplet
                    amountToErode = min((sedimentCapacity - sediment) * erodeSpeed, -deltaHeight)

                    # Use erosion brush to erode from all nodes inside the droplet's erosion radius
                    if erosionRadius <= dropletIndex[1] <= sx - erosionRadius && erosionRadius <= dropletIndex[2] <= sy - erosionRadius

                        brushTranslated = warp(copy(erosionBrush), Translation(-dropletIndex[1],-dropletIndex[2]))
                        #summary(brushTranslated)
                        brushX = axes(brushTranslated, 1)
                        brushY = axes(brushTranslated, 2)
                        #println(brushX)
                        affectedMap = mapEroded[brushX, brushY]
                        #println(affectedMap)
                        afterErosionAffected = affectedMap - (amountToErode * brushTranslated)

                        # clamp erosion to not go below 0
                        afterErosionAffected = max.(0, afterErosionAffected)

                        sediment += sum(affectedMap-afterErosionAffected)
                        mapEroded[brushX, brushY] = afterErosionAffected

                    end
                end

                #println(deltaHeight)
                # Update droplet's speed and water content
                speed = sqrt(max(0,speed^2 - deltaHeight * gravity - friction*speed))
                water *= (1 - evaporateSpeed)
            end
            #dropletPaths[i] = path
    end

    diffMap = map - mapEroded

    kern = KernelFactors.gaussian((blurRadius,blurRadius))
    blurredDiffMap = imfilter(diffMap, kern)

    return map - ( blurFactor * blurredDiffMap + (1-blurFactor) *diffMap)
end


function plotHeightmap(heightmap, sx, sy)
    println("start plotting")
    # Set up sliders to control lighting attributes
    s1, ambient = textslider(0f0:0.01f0:1f0, "ambient", start = 0.55f0)
    s2, diffuse = textslider(0f0:0.025f0:2f0, "diffuse", start = 0.4f0)
    s3, specular = textslider(0f0:0.025f0:2f0, "specular", start = 0.2f0)
    s4, shininess = textslider(2f0.^(2f0:8f0), "shininess", start = 32f0)

    # Set up (r, θ, ϕ) for lightposition
    s5, radius = textslider(2f0.^(0.5f0:0.25f0:20f0), "light pos r", start = 2f0)
    s6, theta = textslider(0:5:180, "light pos theta", start = 30f0)
    s7, phi = textslider(0:5:360, "light pos phi", start = 45f0)

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
        showaxes=false,
        show_axis=false
    )
    #Makie.scatter!(scene2, map(v -> [v], lp), color=:yellow, markersize=1f0)

    # Combine scene
    scene = Makie.Scene(resolution=(1500, 1200), SSAO=ssao_attrib)
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
        1:sx, 1:sy, 1000*heightmap,
        colormap = :Spectral
    )
    scene = vbox(p1)
    #popdisplay(AbstractPlotting.PlotDisplay())
    AbstractPlotting.inline!(false)
    display(scene)
    println("End plot")


    scene

end

function plotMesh(heightmap)

    meshedPts, connectPts = generateMesh(heightmap)

    # Set up sliders to control lighting attributes
    s1, ambient = textslider(0f0:0.01f0:1f0, "ambient", start = 0.55f0)
    s2, diffuse = textslider(0f0:0.025f0:2f0, "diffuse", start = 0.6f0)
    s3, specular = textslider(0f0:0.025f0:2f0, "specular", start = 0.2f0)
    s4, shininess = textslider(2f0.^(2f0:8f0), "shininess", start = 8f0)

    # Set up (r, θ, ϕ) for lightposition
    s5, radius = textslider(2f0.^(0.5f0:0.25f0:20f0), "light pos r", start = 5f0)
    s6, theta = textslider(0:5:180, "light pos theta", start = 55f0)
    s7, phi = textslider(0:5:360, "light pos phi", start = 95f0)

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
    #Makie.scatter!(scene2, map(v -> [v], lp), color=:yellow, markersize=1f0)

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




println("Begin")

sx = 512
sy = 512
hScale = 2.0
seed = 12313
numIterations = 200000
Random.seed!(seed)
erosionBrush = initializeBrush(erosionRadius)
depositionBrush = initializeBrush(erosionRadius)
heightmap = generateNoise(sx, sy, seed)

#dropPaths = Array{Matrix{Float64}}(undef, numIterations)
image2 = erodeTerrain(heightmap, numIterations)
#plot(image)

#r = range(-10, 10, length=1000)
#zs = [sin(2x+y) for x in r, y in r]

plotHeightmap(image2, sx, sy)
#plotMesh(image2)

#gr(legend=false)
#plt = Plots.heatmap(image2)
#plt = Plots.surface(image2, camera=(20,50))
# for path in dropPaths
#     x = path[:, 1]
#     y = path[:, 2]
#
#     plot!(x, y, linecolor=:white)
#     scatter!(x[1:1], y[1:1], markercolor=:white)
# end
#
# rand_cnt = 100
# x = (sx - 1)*rand(rand_cnt) + ones(rand_cnt)
# y = (sy - 1)*rand(rand_cnt) + ones(rand_cnt)
# ux = Array{Float64}(undef, rand_cnt)
# uy = Array{Float64}(undef, rand_cnt)
# gx = Array{Float64}(undef, rand_cnt)
# gy = Array{Float64}(undef, rand_cnt)
#
#
# gx_t, gy_t = GR.gradient(1:sx, 1:sy, heightmap)
#
#
# for i = 1:rand_cnt
#     strc = calcHeightAndGradient(heightmap, x[i], y[i])
#     ux[i] = strc.gradientX
#     uy[i] = strc.gradientY
#
#     x_int::Int = floor(x[i])
#     y_int::Int = floor(y[i])
#     gx[i] = gx_t[x_int]
#     gy[i] = gy_t[y_int]
#
# end
#quiver!(x, y, quiver=(ux, uy))

#display(plt)
println("done")


end
