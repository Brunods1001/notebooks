### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 7f6ae0f5-e682-47b2-ac10-76548fa8330c
using Unitful

# ╔═╡ ac580b2c-fafc-11eb-3e6c-f9b2c8deb462
md"""
# Numerical methods

Solving for a falling body.

"""

# ╔═╡ e3573ca8-87bd-452f-a5b9-a14fadd385a7
md"""
![falling](https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcT3x918lP6IO8R53VUWM4Rcby4yxdxL97ju2Vur8CCcsuBB1zdUzz1HVcBQT1caLOj0-6g&usqp=CAU)
"""

# ╔═╡ b5dd18b4-e016-49f9-a65b-6eee57dc3fa3
md"""
$$\sum{F} = ma$$
$$a = \frac{\sum{F}}{m}$$
$$\frac{dv}{dt} = \frac{\sum{F}}{m}$$
$$\sum{F} = F_U - F_D$$

The force of gravity is:

$$F_D = mg$$
Assume drag is linearly proportional to velocity.

$$F_U = cv$$

Then

$$\frac{dv}{dt} = \frac{cv}{m} - g$$

Where the differential may be approximated as a finite difference:

$$\frac{dv}{dt} \approx \frac{v(t_2) - v(t_1)}{t_2 - t_1} = \frac{cv}{m} - g$$

Which can be rearranged:

$$v(t_2) = v(t_1) + (\frac{cv}{m} - g) * (t_2 - t_1)$$

"""

# ╔═╡ 7d5cb026-c998-4d64-878e-7c15c5ce1130
begin
	m = 68.1u"kg"
	g = 9.81u"m/s^2"
	c = 12.5u"kg/s"
end

# ╔═╡ daa6cc50-315e-450c-8f03-d5e906e610f2
Δt = 2u"s"

# ╔═╡ 1d53a572-c7f7-4115-a5ce-30467387a2d5
v₀ = 0u"m/s"

# ╔═╡ 60835fec-b461-4d80-81dc-88f9006540fa
v₁ = v₀ + (g - c * v₀ / m) * Δt

# ╔═╡ 5c16ece4-b2bf-42f5-b347-bf04bbb04bc9
g, c, m

# ╔═╡ c3f7aeb9-e9e8-4b51-a34e-746113fab553
v₂ = v₁ + (g - c * v₁ / m) * Δt

# ╔═╡ 79eea902-419b-4602-baa1-f7512956bc9e
v₃ = v₂ + (g - c * v₂ / m) * Δt

# ╔═╡ 322b3133-0215-41a0-998c-68e8951a5799
v(t₂, Δt = 1u"s") = t₂ <= 0u"s" ? 0u"m/s" : v(t₂ - Δt, Δt) + (g - c * v(t₂ - Δt, Δt) / m) * (Δt)

# ╔═╡ df2536f6-73f5-4f01-be1d-791863f1e736
v(10u"s", 0.5u"s")

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
Unitful = "~1.9.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a981a8ef8714cba2fd9780b22fd7a469e7aaf56d"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.9.0"
"""

# ╔═╡ Cell order:
# ╟─ac580b2c-fafc-11eb-3e6c-f9b2c8deb462
# ╟─e3573ca8-87bd-452f-a5b9-a14fadd385a7
# ╟─b5dd18b4-e016-49f9-a65b-6eee57dc3fa3
# ╠═7f6ae0f5-e682-47b2-ac10-76548fa8330c
# ╠═7d5cb026-c998-4d64-878e-7c15c5ce1130
# ╠═daa6cc50-315e-450c-8f03-d5e906e610f2
# ╠═1d53a572-c7f7-4115-a5ce-30467387a2d5
# ╠═60835fec-b461-4d80-81dc-88f9006540fa
# ╠═5c16ece4-b2bf-42f5-b347-bf04bbb04bc9
# ╠═c3f7aeb9-e9e8-4b51-a34e-746113fab553
# ╠═79eea902-419b-4602-baa1-f7512956bc9e
# ╠═322b3133-0215-41a0-998c-68e8951a5799
# ╠═df2536f6-73f5-4f01-be1d-791863f1e736
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
