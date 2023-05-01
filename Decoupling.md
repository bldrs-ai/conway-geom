# Decoupling IFC Parsing and Geometry Processing in Web-IFC WASM Module

## Introduction

The Web-IFC WASM module is an essential tool for enabling viewing and editing of the Industry Foundation Classes (IFC) format on the Web. To address the need for fast access to streaming assets and updates, we have decided to make significant modifications to the existing Web-IFC WASM codebase. This document explains the reasons and decisions behind these changes.

## Objectives

The primary objectives of these modifications are:

* Decoupling IFC parsing and geometry processing to create a more modular system.

* Implementing a custom parsing and indexing system for faster access to streaming assets and updates.

* Extending the WASM module to handle geometry extraction and transformations for the IFC format and other 3D file formats.

* Developing a querying system for fast and efficient communication between the system and the WASM module.


## Decoupling IFC Parsing and Geometry Processing

The current Web-IFC WASM module intertwines IFC parsing and geometry processing, which limits flexibility and adaptability to different use cases. By decoupling these processes, we aim to create a more modular and efficient system that can cater to various project requirements.

## Custom Parsing and Indexing System

To enable fast access to streaming assets and updates, we will develop a custom parsing and indexing implementation. This system will allow for efficient data retrieval and updates as the module processes the IFC files. By using this implementation, we can optimize performance and provide a more responsive user experience.

## Extending WASM Module for Geometry Extraction and Transformations

We will extend the WASM module to handle geometry extraction and transformations for the IFC format and other 3D file formats we plan to support. By exposing these functions via Emscripten, the system will be able to perform all necessary transforms and extractions required for different 3D file formats.

## Querying System for Fast Communication

Post indexing, the system will communicate with the WASM module through a querying system. This system will enable users to make fast queries on subsections of the full IFC file. For example, if a user wants the geometry for just the walls or floors of a model, they can send a query to the WASM module, which will return the requested data in a fast and efficient manner.

## Contributing to the IFCjs Team

As part of the project, we are considering contributing their modified Web-IFC WASM codebase, including the parsing and indexing system, to the IFCjs team. This contribution will help improve the IFCjs ecosystem and provide more options for developers and users working with IFC files and other 3D formats.

## Conclusion

The proposed changes to the Web-IFC WASM module aim to create a more modular and efficient system that can handle the requirements of various projects. By decoupling IFC parsing and geometry processing, implementing a custom parsing and indexing system, extending the WASM module, and developing a querying system, the team will provide a more responsive and flexible solution for collaborative viewing and editing of 3D file formats.
