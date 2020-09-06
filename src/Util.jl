module Util

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

export HeightAndGradient

function generateNoise(sx::Int, sy::Int, seed::Int; freq=0.005, d::Int=4, pers=0.3)

        initnoise(seed)


        image = Array{Float32}(undef, sx, sy)

        for pos in Iterators.product(1:sx, 1:sy)
                ns = noise(pos[1]*freq, pos[2]*freq,
                        detail=d,
                        persistence=pers)
                image[pos[1], pos[2]] = ns
        end


        #image_ridge = fill(1.0, size(image)) - 2*abs.(image - fill(0.5, size(image)))
        return image
end

function generateMesh(heightmap, hScale, interpolation=1)
    sz = size(heightmap)

    meshedPts = Array{Float32}(undef, sz[1]*sz[2], 3)

    connSz = 2*(sz[1]-1)*(sz[2]-1)
    connectPts = Array{Int32}(undef, connSz, 3)

    for (i,j) in Iterators.product(1:sz[1], 1:sz[2])
        meshedPts[(j-1)*sz[1] + i, :] = [i  j  hScale*heightmap[i, j]]

        #println(i, j)

        if i < sz[1] && j < sz[2]
            idx = (j-1)*(sz[1]-1) + i
            connectPts[2*(idx - 1) + 1, :] =        [(j-1)*sz[1] + i,     (j-1)*sz[1] + i+ 1,  (j)*sz[1] + i]
            connectPts[2*(idx - 1) + 2, :] =    [(j-1)*sz[1] + i + 1, (j)*sz[1] + i + 1,    (j)*sz[1] + i]

            #println(connectPts)
        end
    end


    return meshedPts, connectPts

end

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

function centralGradient(z)

        arraySize = size(z)

        g = zeros(arraySize[1], arraySize[2], 2)

        # x gradient
        g[1:end-1, :, 1] = z[2:end, :]
        g[2:end, :, 1] -= z[1:end-1, :]

        g[1, :, 1] -= z[1, :]
        g[end, :, 1] += z[end, :]

        g[2:end-1, :, 1] *= 0.5

        # y gradient
        g[:, 1:end-1, 2] = z[:, 2:end]
        g[:, 2:end, 2] -= z[:, 1:end-1]

        g[:, 1, 2] -= z[:, 1]
        g[:, end, 2] += z[:, end]

        g[:, 2:end-1, 2] *= 0.5

        return g

end

function cosLocalSlope(z)

        gradient = centralGradient(z)
        sizeX = size(z,1)
        sizeY = size(z,2)

        cosSlope::Array{Float32} = (sqrt.(gradient[:,:,1].^2 + gradient[:,:,2].^2 + ones(sizeX, sizeY, 1))).^-1

        return dropdims(cosSlope, dims=3)
end

function sinLocalSlope(z)

        gradient = centralGradient(z)
        sizeX = size(z,1)
        sizeY = size(z,2)

        sinSlope::Array{Float32} = sqrt.(ones(sizeX, sizeY, 1) - (gradient[:,:,1].^2 + gradient[:,:,2].^2 + ones(sizeX, sizeY, 1)).^-1)

        return dropdims(sinSlope, dims=3)
end

function normal(z)
        gradient = centralGradient(z)
        result = ones(Float32,size(z,1), size(z,2), 3)
        result[:,:, 1:2] = gradient

        result = result ./ sqrt.(gradient[:,:, 1].^2 .+ gradient[:,:, 2].^2 .+ 1)

        return result
end

function projectionOnPlane(x, n)

A =     [1 - n[1]        -n[1]*n[2]      -n[1]*n[3];
        -n[1]*n[2]      1 - n[2]        -n[2]*n[3];
        -n[1]*n[3]      -n[2]*n[3]      1 - n[3]  ]

return A * x


end



function linearInterpolation(posX, posY, z)
        #float_type = typeof(z[1,1])


        sz = size(posX)
        one = 1

        coordX = floor.(Int, posX)
        coordY = floor.(Int, posY)

        # Calculate droplet's offset inside the cell (0,0) = at NW node, (1,1) = at SE node
        x = posX - coordX
        y = posY - coordY

        # Calculate heights of the four nodes of the droplet's cell
        z_NW = z[coordX, coordY]
        z_NE = z[coordX+one, coordY]
        z_SW = z[coordX, coordY+one]
        z_SE = z[coordX+one, coordY+one]

        # Calculate height with bilinear interpolation of the heights of the nodes of the cell
        z_interp = z_NW .* (one - x) .* (one - y) + z_NE .* x .* (one - y) + z_SW .* (one - x) .* y + z_SE .* x .* y;

        return z_interp
end


# z = rand(4,4)
# x = [1.23 1.25]
# y = [3.25 3.12]
# #interp = linearInterpolation([1.23 1.25], [3.25 3.12], z)
# #println(interp)
# #println(size(interp))
#
# using Interpolations
#
# itp = interpolate(z, BSpline(Linear()))
# result = itp.(x, y)
# println(result)

# z = [
# 1 4 4 1;
# 1 4 4 1;
# 1 4 4 1;
# 1 4 4 1;
# ]
#
# normZ = normal(z)
# display(normZ)
#
# v = Array{Tuple}(undef, size(z))
#
# v[:,:] .= (1,2,3)
#
#
# # result = Array{Float32}(undef,size(z,1),size(z,2),3)
# #
# # for (i,j) in Iterators.product(1:size(z,1), 1:size(z,2))
# #         result[i,j,:] = projectionOnPlane(v[i,j,:], normZ[i,j,:])
# # end
#
# result = projectionOnPlane.(v, normZ)
#
# display(result)
end
