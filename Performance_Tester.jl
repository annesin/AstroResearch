using Pkg
Pkg.clone("https://github.com/JuliaGraphics/Gtk.jl.git")
Pkg.add("Gtk")
using Gtk
Pkg.clone("https://github.com/timholy/ProfileView.jl.git")
Pkg.add("ProfileView")

function myfunc()
    A = rand(200, 200, 400)
    maximum(A)
end

myfunc()

using Profile
Profile.clear()
@profile myfunc()
