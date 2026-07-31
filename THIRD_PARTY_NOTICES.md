# Third-Party Notices

This file records third-party software, data, models, assets, and other materials distributed with MetaSpici or required by a MetaSpici release.

## Current disclosure status

No source tree, dependency manifest, lockfile, container definition, vendored-code inventory, or accessible contents of the target GitHub repository were available in the materials used to prepare this notice. Consequently, no third-party component can be affirmatively identified, versioned, or assigned a license from the supplied materials alone.

This statement is not a representation that MetaSpici has no third-party dependencies. It means that a dependency-level notice cannot responsibly be completed until the exact release contents are audited.

## Confirmed bundled components

No bundled third-party component was identified in the supplied release materials.

Scientific methods and platform names mentioned in the patent materials, including Hi-C, Meta-Hi-C, 3C, Micro-C, Pore-C, SPRITE, and GAM, are not treated as bundled software dependencies merely because they are referenced.

## Required release audit

Before publishing or tagging any MetaSpici release, a maintainer must:

1. inventory direct and transitive runtime, build, test, documentation, workflow, container, and optional dependencies from the exact release commit;
2. inventory vendored source, copied snippets, generated code, models, reference databases, example datasets, fonts, icons, images, and documentation;
3. record each component's canonical name, exact version or commit, source URL, copyright notice, SPDX license identifier, and role in MetaSpici;
4. preserve all license texts, attribution notices, source-offer requirements, modification notices, and other obligations required by those licenses;
5. confirm that every dependency license is compatible with the intended noncommercial and commercial distribution models;
6. separate development-only tools from components actually distributed to users;
7. generate and retain a software bill of materials where practical; and
8. update this file in the same pull request as any dependency change that alters notice obligations.

## Notice format for each component

Add one section per component using this form:

```text
### Component name

- Version or commit:
- Source:
- Copyright:
- License:
- Used for:
- Distributed with MetaSpici: yes/no
- Required notices or modifications:
```

Include the component's complete license text when its terms require reproduction of that text. Do not rely only on a package-registry label; verify the license against the source distribution for the exact version.

## No relicensing of third-party materials

The MetaSpici license applies only to material that the MetaSpici licensor has authority to license. Third-party materials remain governed by their respective licenses and notices. Nothing in `LICENSE`, `COMMERCIAL_LICENSE.md`, or a MetaSpici commercial agreement overrides third-party terms unless the relevant third-party rights holder expressly agrees in writing.

Questions about a third-party notice should be raised through the GitHub issue tracker before redistribution.
