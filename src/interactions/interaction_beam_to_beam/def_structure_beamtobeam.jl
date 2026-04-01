# Structure to store the contact information between two beams in contact
struct BeamPairInContact
    beam1::Int
    beam2::Int
    ξᶜ₁::Float64
    ξᶜ₂::Float64
    xᵖ₁::Vec3
    xᵖ₂::Vec3
    G1P::Mat312
    G2P::Mat312
    gᵀ₁::Float64
    gᵀ₂::Float64
    gᵀᴾ₁::Float64
    gᵀᴾ₂::Float64
end

# Container structure to store all beam-to-beam contacts (mutable to allow updates)
mutable struct Beam2BeamContacts
    contacts::Dict{Tuple{Int,Int}, BeamPairInContact}
end

# Constructor to initialize an empty contact container
function Beam2BeamContacts()
    return Beam2BeamContacts(Dict{Tuple{Int,Int}, BeamPairInContact}())
end

# Function to add or update a contact in the container
function Update_BeamPair!(all_contacts::Beam2BeamContacts, contact::BeamPairInContact)
    key = (contact.beam1, contact.beam2)
    all_contacts.contacts[key] = contact
end

# Function to check if a contact exists
function Check_BeamPair(all_contacts::Beam2BeamContacts, beam1_ind::Int, beam2_ind::Int)
    key = (beam1_ind, beam2_ind)
    return haskey(all_contacts.contacts, key)
end

# Function to get a contact
function Get_BeamPair(all_contacts::Beam2BeamContacts, beam1_ind::Int, beam2_ind::Int)
    key = (beam1_ind, beam2_ind)
    return all_contacts.contacts[key]
end

# Function to clear all contacts from the container
function Clear_All_Contacts!(all_contacts::Beam2BeamContacts)
    empty!(all_contacts.contacts)
end

