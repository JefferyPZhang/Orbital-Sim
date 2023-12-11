# CubeSat Orbital Simulator

Welcome to this very crude CubeSat orbital simulator! Here are a couple of things to keep in mind as you get the simulation running:
Position is defined using spherical coordinates, where  is the distance from the center of the Earth, phi is the inclination angle (fully-defined from 0 to pi) and θ is the azimuthal angle (fully-defined from 0 to 2pi).
Since we only want to consider circular orbits, the velocity vector magnitude will always be kept at: v = sqrt(2GM/r)

You will be able to set the mass of the satellite (1U, 2U, 3U, 6U, 12U), where each unit of the satellite is roughly 1 kg.
You will be able to set the initial position of the satellite (in spherical coordinates - for rho, you only need to specify the altitude at which the satellite is orbiting, which must be between 150 km and 2000 km above the surface) and the directional vector of the velocity, with the "zero" in reference to the top of the satellite facing an extnesion of the north pole (while keeping tangent to the Earth's surface).
As the simulation is running, you will be able to use specify a new orbital altitude and the satellite will immediately adjust using Hohmann transfers (again, must be between 150 km and 2000 km above the Earth's surface). You can also adjust the playback speed of the simulation.
