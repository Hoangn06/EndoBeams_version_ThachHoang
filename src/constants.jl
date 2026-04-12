const Mat22{Float64} = SMatrix{2,2,Float64,4}

const Mat33{Float64} = SMatrix{3,3,Float64,9}
const Mat34{Float64} = SMatrix{3,4,Float64,12}
const Mat37{Float64} = SMatrix{3,7,Float64,21}
const Mat67{Float64} = SMatrix{6,7,Float64,42}
const Mat66{Float64} = SMatrix{6,6,Float64,36}
const Mat77{Float64} = SMatrix{7,7,Float64,49}

const Mat312{Float64} = SMatrix{3,12,Float64,36}
const Mat224{Float64} = SMatrix{2,24,Float64,48}
const Mat324{Float64} = SMatrix{3,24,Float64,72}


const Vec1{Float64} = SVector{1,Float64}
const Vec2{Float64} = SVector{2,Float64}
const Vec3{Float64} = SVector{3,Float64}
const Vec4{Float64} = SVector{4,Float64}
const Vec5{Float64} = SVector{5,Float64}
const Vec6{Float64} = SVector{6,Float64}
const Vec7{Float64} = SVector{7,Float64}
const Vec10{Float64} = SVector{10,Float64}

const Vec24{Float64} = SVector{24,Float64} 

const ID3 = Diagonal(Vec3(1,1,1))
const ID6 = Diagonal(Vec6(1.0,1.0,1.0,1.0,1.0,1.0))
const Idev = Diagonal(Vec3(2/3, 1.0, 1.0))

const I = Vec3(1.0, 0.0, 0.0)

